find . -maxdepth 1 -name ".svn" -prune -o -wholename "." -o -type d -print -exec touch '{}'/__init__.py \;
