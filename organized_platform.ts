const { promises: { copyFile, readFile, stat, access}, constants } = require('fs');
const {join} = require('path');


const metadataSourcePath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\metadata.cart.2022-06-08.json";
const sourceFolderPath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\dna methylation"


type MetadataRecord = {
    file_name: string;
    platform: string;
};

const dataPlatformDestination: { [key: string]: string; } = {
    "Illumina Human Methylation 27": "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\dna methylation by platform\\HM27",
    "Illumina Human Methylation 450": "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\dna methylation by platform\\HM450",

};

(async () => {
    //todo:read the file in as a string
const metadataFileContents = await readFile(metadataSourcePath, 'utf8');
    
//todo: pares the string into a JSON object
    const metadata: MetadataRecord[] = JSON.parse(metadataFileContents);
    
    //todo: write the JSON object to the console
    //console.log(metadata);

    for (const record of metadata){
    try{
        const {file_name, platform } = record;

        
        
        const sourceFilePath: string = join(sourceFolderPath, file_name);

        await access(sourceFilePath, constants.R_OK);
        
        const destinationFolderPath: string = dataPlatformDestination[platform];
        
        if(!destinationFolderPath) {
            throw new Error (`No destination folder for ${platform}`);
        }
        //todo: verify source file exist       
           
        //todo: get full destination path
        const destinationFilePath: string = join(destinationFolderPath, file_name);

        //todo: verify destination folder exists
            
        //todo: verify write permission to destination folder
       
    

        //todo: copy source file to new destination
        await access(destinationFolderPath, constants.W_OK);

        await copyFile(sourceFilePath, destinationFilePath);
        console.log("copyFile")
    } catch(error){
        console.log(error)
        console.log(`Failed to process ${record.platform}`)
    }

        

        
    }

})();