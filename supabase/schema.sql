-- Thorntales: profiles, comments, avatars.
-- Run once in Supabase Dashboard -> SQL Editor -> New query -> Run.
-- Safe to re-run: everything is "create if not exists" / "or replace" / drop-and-recreate for policies.

-- =========================================================
-- Profiles: one row per user (name, avatar, admin flag)
-- =========================================================

create table if not exists public.profiles (
  id           uuid primary key references auth.users (id) on delete cascade,
  display_name text not null check (char_length(display_name) between 1 and 40),
  avatar_url   text check (avatar_url is null or char_length(avatar_url) <= 500),
  is_admin     boolean not null default false,
  created_at   timestamptz not null default now()
);

alter table public.profiles enable row level security;

-- Create a profile automatically when a user is added in Authentication -> Users.
-- Default name = the part of the email before "@"; the user can change it on the account page.
create or replace function public.handle_new_user()
returns trigger
language plpgsql
security definer
set search_path = ''
as $$
begin
  insert into public.profiles (id, display_name)
  values (new.id, left(split_part(new.email, '@', 1), 40))
  on conflict (id) do nothing;
  return new;
end;
$$;

drop trigger if exists on_auth_user_created on auth.users;
create trigger on_auth_user_created
  after insert on auth.users
  for each row execute function public.handle_new_user();

-- Profiles for users that already exist
insert into public.profiles (id, display_name)
select id, left(split_part(email, '@', 1), 40) from auth.users
on conflict (id) do nothing;

-- Helper: is the current user an admin? (security definer avoids RLS recursion)
create or replace function public.is_admin()
returns boolean
language sql
stable
security definer
set search_path = ''
as $$
  select coalesce((select is_admin from public.profiles where id = auth.uid()), false);
$$;

-- Everyone can read profiles (names/avatars are shown next to comments)
drop policy if exists "profiles are public" on public.profiles;
create policy "profiles are public" on public.profiles
  for select using (true);

-- A user can edit only their own profile...
drop policy if exists "users edit own profile" on public.profiles;
create policy "users edit own profile" on public.profiles
  for update to authenticated
  using (id = auth.uid())
  with check (id = auth.uid());

-- ...and only these columns (is_admin can be changed only here, in the SQL editor)
revoke insert, update, delete on public.profiles from anon, authenticated;
grant select on public.profiles to anon, authenticated;
grant update (display_name, avatar_url) on public.profiles to authenticated;


-- =========================================================
-- Comments
-- thread = what the comment belongs to, e.g. 'poem/salvia/2018-06-24-обойди-меня-за-девять-верст'
-- =========================================================

create table if not exists public.comments (
  id         bigint generated always as identity primary key,
  thread     text not null check (char_length(thread) between 1 and 300),
  user_id    uuid not null default auth.uid() references public.profiles (id) on delete cascade,
  body       text not null check (char_length(btrim(body)) between 1 and 4000),
  created_at timestamptz not null default now(),
  updated_at timestamptz
);

create index if not exists comments_thread_created_idx on public.comments (thread, created_at);

alter table public.comments enable row level security;

create or replace function public.touch_comment()
returns trigger
language plpgsql
set search_path = ''
as $$
begin
  new.updated_at := now();
  return new;
end;
$$;

drop trigger if exists comments_touch on public.comments;
create trigger comments_touch
  before update on public.comments
  for each row execute function public.touch_comment();

-- Everyone (also not logged in) can read comments
drop policy if exists "comments are public" on public.comments;
create policy "comments are public" on public.comments
  for select using (true);

-- Logged-in users write comments only as themselves
drop policy if exists "users add own comments" on public.comments;
create policy "users add own comments" on public.comments
  for insert to authenticated
  with check (user_id = auth.uid());

-- Only the author edits a comment
drop policy if exists "authors edit own comments" on public.comments;
create policy "authors edit own comments" on public.comments
  for update to authenticated
  using (user_id = auth.uid())
  with check (user_id = auth.uid());

-- The author or an admin deletes a comment
drop policy if exists "authors or admin delete comments" on public.comments;
create policy "authors or admin delete comments" on public.comments
  for delete to authenticated
  using (user_id = auth.uid() or public.is_admin());

revoke insert, update, delete on public.comments from anon, authenticated;
grant select on public.comments to anon, authenticated;
grant insert (thread, body) on public.comments to authenticated;
grant update (body) on public.comments to authenticated;
grant delete on public.comments to authenticated;


-- =========================================================
-- Avatars: public bucket, each user writes only into their own folder <user id>/...
-- =========================================================

insert into storage.buckets (id, name, public, file_size_limit, allowed_mime_types)
values ('avatars', 'avatars', true, 1048576, array['image/jpeg', 'image/png', 'image/webp'])
on conflict (id) do update
  set public = excluded.public,
      file_size_limit = excluded.file_size_limit,
      allowed_mime_types = excluded.allowed_mime_types;

drop policy if exists "avatars upload own" on storage.objects;
create policy "avatars upload own" on storage.objects
  for insert to authenticated
  with check (bucket_id = 'avatars' and (storage.foldername(name))[1] = auth.uid()::text);

drop policy if exists "avatars update own" on storage.objects;
create policy "avatars update own" on storage.objects
  for update to authenticated
  using (bucket_id = 'avatars' and (storage.foldername(name))[1] = auth.uid()::text);

drop policy if exists "avatars delete own" on storage.objects;
create policy "avatars delete own" on storage.objects
  for delete to authenticated
  using (bucket_id = 'avatars' and (storage.foldername(name))[1] = auth.uid()::text);

-- Listing own folder is needed to remove old avatar files
drop policy if exists "avatars list own" on storage.objects;
create policy "avatars list own" on storage.objects
  for select to authenticated
  using (bucket_id = 'avatars' and (storage.foldername(name))[1] = auth.uid()::text);


-- =========================================================
-- Poem edits made on the site by the author (admin).
-- One row = the whole poem as it should be shown now:
--   poem_id found in data/poetry/<catalog>.txt -> replaces that poem;
--   poem_id not in the file                    -> a new poem in that catalog.
-- Deleting a row brings back the text from the file.
-- =========================================================

create table if not exists public.poem_edits (
  catalog    text not null check (char_length(catalog) between 1 and 80),
  poem_id    text not null check (char_length(poem_id) between 1 and 100),
  title      text not null check (char_length(title) between 1 and 300),
  date       text not null default '' check (char_length(date) <= 40),
  body       text not null check (char_length(body) <= 50000),
  tags       text[] not null default '{}' check (cardinality(tags) <= 30),
  updated_at timestamptz not null default now(),
  primary key (catalog, poem_id)
);

alter table public.poem_edits enable row level security;

drop trigger if exists poem_edits_touch on public.poem_edits;
create trigger poem_edits_touch
  before update on public.poem_edits
  for each row execute function public.touch_comment();

drop policy if exists "poem edits are public" on public.poem_edits;
create policy "poem edits are public" on public.poem_edits
  for select using (true);

drop policy if exists "admin writes poem edits" on public.poem_edits;
create policy "admin writes poem edits" on public.poem_edits
  for all to authenticated
  using (public.is_admin())
  with check (public.is_admin());

revoke insert, update, delete on public.poem_edits from anon, authenticated;
grant select on public.poem_edits to anon, authenticated;
grant insert, update, delete on public.poem_edits to authenticated;


-- =========================================================
-- Collection edits made on the site by the admin: name, description, small picture.
-- A row overrides the entry with the same id in data/poetry/catalogs.json.
-- =========================================================

create table if not exists public.catalog_edits (
  catalog    text primary key check (char_length(catalog) between 1 and 80),
  label      text not null check (char_length(label) between 1 and 200),
  descr      text not null default '' check (char_length(descr) <= 4000),
  image_url  text check (image_url is null or char_length(image_url) <= 500),
  updated_at timestamptz not null default now()
);

alter table public.catalog_edits enable row level security;

drop trigger if exists catalog_edits_touch on public.catalog_edits;
create trigger catalog_edits_touch
  before update on public.catalog_edits
  for each row execute function public.touch_comment();

drop policy if exists "catalog edits are public" on public.catalog_edits;
create policy "catalog edits are public" on public.catalog_edits
  for select using (true);

drop policy if exists "admin writes catalog edits" on public.catalog_edits;
create policy "admin writes catalog edits" on public.catalog_edits
  for all to authenticated
  using (public.is_admin())
  with check (public.is_admin());

revoke insert, update, delete on public.catalog_edits from anon, authenticated;
grant select on public.catalog_edits to anon, authenticated;
grant insert, update, delete on public.catalog_edits to authenticated;

-- Pictures of collections: public bucket, only the admin uploads
insert into storage.buckets (id, name, public, file_size_limit, allowed_mime_types)
values ('catalog-images', 'catalog-images', true, 1048576, array['image/jpeg', 'image/png', 'image/webp'])
on conflict (id) do update
  set public = excluded.public,
      file_size_limit = excluded.file_size_limit,
      allowed_mime_types = excluded.allowed_mime_types;

drop policy if exists "catalog images admin insert" on storage.objects;
create policy "catalog images admin insert" on storage.objects
  for insert to authenticated
  with check (bucket_id = 'catalog-images' and public.is_admin());

drop policy if exists "catalog images admin update" on storage.objects;
create policy "catalog images admin update" on storage.objects
  for update to authenticated
  using (bucket_id = 'catalog-images' and public.is_admin());

drop policy if exists "catalog images admin delete" on storage.objects;
create policy "catalog images admin delete" on storage.objects
  for delete to authenticated
  using (bucket_id = 'catalog-images' and public.is_admin());

drop policy if exists "catalog images admin list" on storage.objects;
create policy "catalog images admin list" on storage.objects
  for select to authenticated
  using (bucket_id = 'catalog-images' and public.is_admin());


-- =========================================================
-- Who sees poems with a tag.
--   no row / mode 'public' -> everyone
--   mode 'signed_in'       -> everyone who is logged in
--   mode 'listed'          -> only users in tag_readers (empty list = only the admin)
-- The admin sees everything. A poem with several restricted tags is shown only
-- if every one of them lets the reader in.
-- NOTE: this protects poems stored in poem_edits. A poem that is also in a txt file
-- is only hidden on the pages: the file itself is public on GitHub.
-- =========================================================

create table if not exists public.tag_access (
  tag        text primary key check (char_length(tag) between 1 and 40),
  mode       text not null default 'public' check (mode in ('public', 'signed_in', 'listed')),
  updated_at timestamptz not null default now()
);

create table if not exists public.tag_readers (
  tag     text not null references public.tag_access (tag) on delete cascade on update cascade,
  user_id uuid not null references public.profiles (id) on delete cascade,
  primary key (tag, user_id)
);

alter table public.tag_access enable row level security;
alter table public.tag_readers enable row level security;

drop trigger if exists tag_access_touch on public.tag_access;
create trigger tag_access_touch
  before update on public.tag_access
  for each row execute function public.touch_comment();

drop policy if exists "admin manages tag access" on public.tag_access;
create policy "admin manages tag access" on public.tag_access
  for all to authenticated
  using (public.is_admin())
  with check (public.is_admin());

drop policy if exists "admin manages tag readers" on public.tag_readers;
create policy "admin manages tag readers" on public.tag_readers
  for all to authenticated
  using (public.is_admin())
  with check (public.is_admin());

revoke all on public.tag_access, public.tag_readers from anon, authenticated;
grant select, insert, update, delete on public.tag_access, public.tag_readers to authenticated;

-- Can the current reader see a poem with these tags?
create or replace function public.can_see_tags(tags text[])
returns boolean
language sql
stable
security definer
set search_path = ''
as $$
  select public.is_admin() or not exists (
    select 1
    from public.tag_access a
    where a.tag = any (tags)
      and a.mode <> 'public'
      and not (a.mode = 'signed_in' and auth.uid() is not null)
      and not (a.mode = 'listed' and exists (
        select 1 from public.tag_readers r where r.tag = a.tag and r.user_id = auth.uid()))
  );
$$;

-- Poems the current reader must not see (so the pages also hide their txt versions)
create or replace function public.hidden_poems()
returns table (catalog text, poem_id text)
language sql
stable
security definer
set search_path = ''
as $$
  select e.catalog, e.poem_id
  from public.poem_edits e
  where not public.can_see_tags(e.tags);
$$;

grant execute on function public.can_see_tags(text[]) to anon, authenticated;
grant execute on function public.hidden_poems() to anon, authenticated;

-- Restricted poems are not readable at all for those who may not see them
drop policy if exists "poem edits are public" on public.poem_edits;
drop policy if exists "poem edits by tag access" on public.poem_edits;
create policy "poem edits by tag access" on public.poem_edits
  for select using (public.can_see_tags(tags));


-- =========================================================
-- Texts (stories, tales, chatter, ...) edited or added on the site by the admin.
-- One row = the whole text as it should be shown now, with its translations:
--   text_id found in data/texts/<section>.json -> replaces that text;
--   text_id not in the file                    -> a new text in that section;
--   deleted = true                             -> the text is hidden.
-- tr = { "en": { "title": "...", "body": "<p>...</p>" }, "de": { ... } }
-- Deleting a row brings back the text from the file.
-- =========================================================

create table if not exists public.text_edits (
  section    text not null check (char_length(section) between 1 and 80),
  text_id    text not null check (char_length(text_id) between 1 and 100),
  title      text not null check (char_length(title) between 1 and 300),
  date       text not null default '' check (char_length(date) <= 40),
  lang       text not null default 'ru' check (lang in ('ru', 'en', 'de')),
  body       text not null default '' check (char_length(body) <= 400000),
  tr         jsonb not null default '{}'::jsonb,
  deleted    boolean not null default false,
  updated_at timestamptz not null default now(),
  primary key (section, text_id)
);

alter table public.text_edits enable row level security;

drop trigger if exists text_edits_touch on public.text_edits;
create trigger text_edits_touch
  before update on public.text_edits
  for each row execute function public.touch_comment();

drop policy if exists "text edits are public" on public.text_edits;
create policy "text edits are public" on public.text_edits
  for select using (true);

drop policy if exists "admin writes text edits" on public.text_edits;
create policy "admin writes text edits" on public.text_edits
  for all to authenticated
  using (public.is_admin())
  with check (public.is_admin());

revoke insert, update, delete on public.text_edits from anon, authenticated;
grant select on public.text_edits to anon, authenticated;
grant insert, update, delete on public.text_edits to authenticated;


-- =========================================================
-- Make yourself admin (can delete any comment).
-- Put your login email below and run just this statement.
-- =========================================================
-- update public.profiles set is_admin = true
-- where id = (select id from auth.users where email = 'YOUR_EMAIL_HERE');
