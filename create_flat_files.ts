const { promises: { readdir, copyFile } } = require('fs');
const {join} = require('path')

const sourcePath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\target aml data";
const destinationPath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\aml data flat";

(async () => {
  const sourceDirectories = await readdir(sourcePath, {withFileTypes: true});

  for(const sourceDirectory of sourceDirectories) {
    const fullSourcePath = join(sourcePath, sourceDirectory.name);
    
    const directoryContents = await readdir(fullSourcePath, {
      withFileTypes: true,
    });

    for (const file of directoryContents){
      const fullFilePath = join(fullSourcePath, file.name);
      
      if (file.isFile()) {
       
        const destinationFilePath = join(destinationPath, file.name);
        
        await copyFile(fullFilePath, destinationFilePath);

      }
        if(file.isDirectory()){
          const nestedDirectoryContents = await readdir(fullFilePath, {
            withFileTypes: true,
          });
      
          for (const nestedFile of nestedDirectoryContents){
            const nestedFullFilePath = join(fullFilePath, nestedFile.name);
            
            if (nestedFile.isFile()) {
             
              const nestedDestinationFilePath = join(destinationPath, nestedFile.name);
              
              await copyFile(nestedFullFilePath, nestedDestinationFilePath);
      
            }
              
            
          }
        }
      
    }
    }
})();