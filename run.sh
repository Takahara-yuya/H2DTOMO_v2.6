numtask=20
folder_path="Result/tmp"
folder_path2="Result"
exe="H2DTOMO"
#delete tmp
find "$folder_path" -mindepth 1 -delete
find "$folder_path" -type f -delete
echo "Folder deleted successfully."
#
mpirun -np $numtask ./$exe
