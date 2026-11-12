#!/bin/bash

set -e

version=$(cat ../VERSION)
revision="unavailable (source archive)"
if [[ -e ../.git ]]; then
    revision=$(git -C .. rev-parse --verify HEAD)
fi

cat > gitver.f90 <<EOF
subroutine gitversion(version,sha)
  implicit none
  character(len=*), intent(out) :: version,sha

  version = '$version'
  sha = '$revision'
end subroutine gitversion
EOF
