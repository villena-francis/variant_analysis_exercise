cd part_2

if md5sum --status -c checksum.txt
then
    echo "All MD5s match."
else
    echo "At least one MD5 does not match."
fi