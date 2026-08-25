set -euo pipefail
ROOT=/root/autodl-tmp/o6u_md_release_3x500ns_v4
echo '--- release top-level archives ---'
find "$ROOT" -maxdepth 2 -type f \( -name '*.tgz' -o -name '*.tar.gz' -o -name '*.zip' \) -printf '%p %s bytes\n' | sort
echo '--- exact source basenames in release ---'
find "$ROOT" -type f \( \
  -name '8KCT_protonated.pdb' -o \
  -name 'O6U_CCD_CGENFF_ATOM_CORRESPONDENCE.tsv' -o \
  -name '8KCT_O6U.xml' -o \
  -name '8KCT_O6U_native_contacts.interactions.normalized.json' -o \
  -name 'step5_input.pdb' -o \
  -name 'minimized.gro' \
\) -printf '%p %s bytes\n' | sort
echo '--- root listing ---'
find "$ROOT" -maxdepth 1 -printf '%f %y %s bytes\n' | sort
