for var in {0..398}
do
    ./build/GoL2L ./data//office0/model_$var.txt ./data/office0/data_$var.txt ./config_line.txt ./data/office0_delta.txt 
done