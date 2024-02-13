for dir in alignment calling genome intervals raw_data
do
    find "$dir" -mindepth 1 -not -name ".gitkeep" -delete
    echo "$dir was emptied"
done
