#!/bin/bash
# Downloads the external inputs into deps/ and instances/:
#   deps/pgl                      Pangolin plane geometry library (header-only), at the commit used for the paper
#   instances/cgshop2019/...      CG:SHOP 2019 challenge instances (the image-based ones are used)
# CPLEX 22.2 must be installed separately (free for academic use through the IBM Academic Initiative).
set -e
cd "$(dirname "$0")"
PGL_COMMIT=d90896ec62432fb73baa33408fa7d5c655cee537
CG_URL=https://cgshop.ibr.cs.tu-bs.de/static/competitions/cgshop_2019/challenge_instances.zip
CG_SHA256=82e0fc3535966ac1c5530cb12b1935b2ffa26f3e04286db1ee8f6487e0488eb0
mkdir -p deps instances
if [ ! -d deps/pgl ]; then
  git clone https://github.com/gfonsecabr/pgl deps/pgl
  git -C deps/pgl checkout -q $PGL_COMMIT
fi
if [ ! -d instances/cgshop2019/challenge_instances ]; then
  curl -fL -o instances/challenge_instances.zip $CG_URL
  echo "$CG_SHA256  instances/challenge_instances.zip" | sha256sum -c -
  mkdir -p instances/cgshop2019
  unzip -q instances/challenge_instances.zip -d instances/cgshop2019
  rm instances/challenge_instances.zip
fi
echo "pgl: $(git -C deps/pgl rev-parse HEAD)"
echo "CG:SHOP 2019: $(ls instances/cgshop2019/challenge_instances/data/images | wc -l) image-based instances"
