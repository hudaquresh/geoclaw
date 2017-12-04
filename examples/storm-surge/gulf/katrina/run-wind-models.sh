for d in */ ; do
    echo "$d"
    cd "$d" 
    make all 
    cd .. 
done


