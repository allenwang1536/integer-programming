#!/usr/bin/env bash
set -euo pipefail

ORTOOLS_VERSION="9.15"
ORTOOLS_BUILD="6755"
ORTOOLS_FULL="${ORTOOLS_VERSION}.${ORTOOLS_BUILD}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENDOR_DIR="$ROOT/.vendor/or-tools/$ORTOOLS_FULL"
TMP_DIR="$ROOT/.vendor/tmp"
mkdir -p "$TMP_DIR"

UNAME_S="$(uname -s)"
ARCH="$(uname -m)"

case "$UNAME_S" in
  Linux)
    if [[ "$ARCH" != "x86_64" ]]; then
      echo "Linux support in this script is currently pinned to x86_64 only."
      exit 1
    fi

    if [[ ! -f /etc/os-release ]]; then
      echo "Cannot detect Linux distro: /etc/os-release missing."
      exit 1
    fi

    # shellcheck disable=SC1091
    source /etc/os-release
    OS_ID="${ID:-}"
    OS_VERSION_ID="${VERSION_ID:-}"

    case "${OS_ID}:${OS_VERSION_ID}" in
      ubuntu:24.04) ASSET="or-tools_amd64_ubuntu-24.04_cpp_v${ORTOOLS_FULL}.tar.gz" ;;
      ubuntu:22.04) ASSET="or-tools_amd64_ubuntu-22.04_cpp_v${ORTOOLS_FULL}.tar.gz" ;;
      ubuntu:20.04) ASSET="or-tools_amd64_ubuntu-20.04_cpp_v${ORTOOLS_FULL}.tar.gz" ;;
      debian:12)    ASSET="or-tools_amd64_debian-12_cpp_v${ORTOOLS_FULL}.tar.gz" ;;
      debian:11)    ASSET="or-tools_amd64_debian-11_cpp_v${ORTOOLS_FULL}.tar.gz" ;;
      *)
        echo "Unsupported Linux distro for this pinned script: ${OS_ID} ${OS_VERSION_ID}"
        exit 1
        ;;
    esac
    ;;

  Darwin)
    # Optional sanity check for local Mac setup
    if ! xcode-select -p >/dev/null 2>&1; then
      echo "Xcode Command Line Tools not found. Run: xcode-select --install"
      exit 1
    fi

    case "$ARCH" in
      x86_64)
        ASSET="or-tools_x86_64_macOS-26.2_cpp_v${ORTOOLS_FULL}.tar.gz"
        ;;
      arm64)
        ASSET="or-tools_arm64_macOS-26.2_cpp_v${ORTOOLS_FULL}.tar.gz"
        ;;
      *)
        echo "Unsupported macOS architecture: $ARCH"
        exit 1
        ;;
    esac
    ;;

  *)
    echo "Unsupported OS: $UNAME_S"
    exit 1
    ;;
esac

URL="https://github.com/google/or-tools/releases/download/v${ORTOOLS_VERSION}/${ASSET}"
ARCHIVE="$TMP_DIR/$ASSET"

if [[ ! -f "$ARCHIVE" ]]; then
  echo "Downloading $URL"
  curl -L --fail --retry 3 -o "$ARCHIVE" "$URL"
fi

if [[ ! -d "$VENDOR_DIR/include" ]]; then
  rm -rf "$VENDOR_DIR"
  mkdir -p "$VENDOR_DIR"
  tar -xzf "$ARCHIVE" -C "$VENDOR_DIR" --strip-components=1
fi

echo "OR-Tools installed at: $VENDOR_DIR"
