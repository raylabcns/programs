% Deletes a directory if its present
function deleteDirectory(dirName)

if exist(dirName,'dir')         
    rmdir(dirName,'s');
    disp(['deleting ' dirName]);
end
