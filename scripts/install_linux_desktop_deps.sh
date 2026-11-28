#!/usr/bin/env bash
# Hosted Ubuntu runner setup only; never called by the Collagene application.
set -euo pipefail

case "${1:-build}" in
  build) collagene_extra_packages=() ;;
  package) collagene_extra_packages=(patchelf xdg-utils) ;;
  *) echo 'Expected build or package dependency profile.' >&2; exit 2 ;;
esac

if [[ -f /etc/apt/sources.list.d/ubuntu.sources ]]; then
  collagene_apt_source=/etc/apt/sources.list.d/ubuntu.sources
elif [[ -f /etc/apt/sources.list ]]; then
  collagene_apt_source=/etc/apt/sources.list
else
  echo 'The Ubuntu runner source list is unavailable.' >&2
  exit 1
fi

# Hosted runners also configure unrelated third-party feeds. Their availability
# must not prevent resolution of these Ubuntu WebView/build packages. Retain
# normal Ubuntu signature and TLS validation; leave the runner configuration alone.
collagene_apt_parts=$(mktemp -d /tmp/collagene-apt-parts.XXXXXX)
trap 'rmdir -- "$collagene_apt_parts"' EXIT
collagene_apt_options=(
  -o "Dir::Etc::sourcelist=$collagene_apt_source"
  -o "Dir::Etc::sourceparts=$collagene_apt_parts"
)
sudo apt-get "${collagene_apt_options[@]}" update
sudo apt-get "${collagene_apt_options[@]}" install --no-install-recommends -y \
  libwebkit2gtk-4.1-dev build-essential curl wget file libxdo-dev libssl-dev \
  libayatana-appindicator3-dev librsvg2-dev pkg-config \
  "${collagene_extra_packages[@]}"
