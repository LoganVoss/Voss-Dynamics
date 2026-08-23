#!/bin/sh

set -eu

script_dir=$(CDPATH= cd -P "$(dirname "$0")" && pwd)
raw_dir="$script_dir/raw"
mkdir -p "$raw_dir"

sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        printf '%s\n' 'error: sha256sum or shasum is required' >&2
        return 127
    fi
}

fetch_to() {
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 --output "$2" "$1"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "$2" "$1"
    else
        printf '%s\n' 'error: curl or wget is required' >&2
        return 127
    fi
}

download_one() {
    name=$1
    url=$2
    expected=$3
    target="$raw_dir/$name"

    if [ -f "$target" ]; then
        actual=$(sha256_of "$target")
        if [ "$actual" = "$expected" ]; then
            printf 'verified existing %s\n' "$name"
            return 0
        fi
        printf 'error: refusing to overwrite %s; expected %s, found %s\n' \
            "$target" "$expected" "$actual" >&2
        return 1
    fi

    temp="$target.part.$$"
    if [ -e "$temp" ]; then
        printf 'error: temporary path already exists: %s\n' "$temp" >&2
        return 1
    fi
    trap 'rm -f "$temp"' 0 1 2 15

    printf 'downloading %s\n' "$name"
    fetch_to "$url" "$temp"
    actual=$(sha256_of "$temp")
    if [ "$actual" != "$expected" ]; then
        printf 'error: checksum mismatch for %s; expected %s, found %s\n' \
            "$name" "$expected" "$actual" >&2
        return 1
    fi

    mv "$temp" "$target"
    trap - 0 1 2 15
    printf 'installed %s\n' "$target"
}

download_one \
    'ECGFiveDays.zip' \
    'https://www.timeseriesclassification.com/aeon-toolkit/ECGFiveDays.zip' \
    '11457a2d590711598eac1a0ab5c58d43c5e3c5c2c86521809d69f7c0b6b3edd1'

download_one \
    'Earthquakes.zip' \
    'https://www.timeseriesclassification.com/aeon-toolkit/Earthquakes.zip' \
    '927cfb732988055850a74efa169dfd633bc7f578095c35319bebb83453501bf1'

printf '%s\n' 'all UCR archives verified'
