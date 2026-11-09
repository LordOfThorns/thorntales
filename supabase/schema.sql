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
-- Make yourself admin (can delete any comment).
-- Put your login email below and run just this statement.
-- =========================================================
-- update public.profiles set is_admin = true
-- where id = (select id from auth.users where email = 'YOUR_EMAIL_HERE');
