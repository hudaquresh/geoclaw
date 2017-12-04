for d in */ ; do
    echo "$d"
    cd "$d"
    source run-wind-models.sh 
    cd ..
done


