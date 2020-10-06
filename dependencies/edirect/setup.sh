#!/bin/bash -norc

DIR="$( cd "$( dirname "$0" )" && pwd )"

cat <<EOF

Trying to establish local installations of any missing Perl modules
(as logged in $DIR/setup-deps.log).
Please be patient, as this step may take a little while.
EOF
mkdir -p "$DIR/_cpan/CPAN"
echo '1;' >> "$DIR/_cpan/CPAN/MyConfig.pm"
perl -I"$DIR/_cpan" "$DIR/setup-deps.pl" </dev/null >"$DIR/setup-deps.log" 2>&1
rm -rf "$DIR/_cpan"

target=bash_profile
if ! grep "$target" "$HOME/.bashrc" >/dev/null 2>&1
then
  if grep 'bashrc' "$HOME/.$target" >/dev/null 2>&1
  then
    target=bashrc
  else
    echo 'source ~/.bash_profile' >> $HOME/.bashrc
  fi
fi
if ! grep "PATH.*edirect" "$HOME/.$target" >/dev/null 2>&1
then
  echo "export PATH=\$PATH:$DIR" >> $HOME/.$target
fi

cd "$DIR"
osname=`uname -s`
cputype=`uname -m`
case "$osname-$cputype" in
  Linux-x86_64 | Darwin-x86_64 )
    ./ftp-cp ftp.ncbi.nlm.nih.gov /entrez/entrezdirect xtract."$osname".gz
    gunzip -f xtract."$osname".gz
  ;;
esac
if [ -f xtract."$osname" ]
then
  chmod +x xtract."$osname"
else
  echo "Unable to download a prebuilt xtract executable; attempting to"
  echo "build one from xtract.go.  A Perl fallback is also available, and"
  echo "will be used if necessary, so please disregard any errors below."
  go build -o xtract."$osname" xtract.go
fi

echo ""
echo "ENTREZ DIRECT HAS BEEN SUCCESSFULLY INSTALLED AND CONFIGURED"
echo ""
