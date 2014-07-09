% checks if a file exist and deltes it.
function deleteFile(fileName)

if exist(fileName,'file')         
    delete(fileName);         
end