const { promises: { copyFile, readFile, stat, access}, constants } = require('fs');
const {join} = require('path');


const metadataSourcePath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\metadata.cart.2022-06-08.json";
const sourceFolderPath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\aml data flat"


type MetadataRecord = {
    file_name: string;
    data_category: string;
};

const dataCategoryDestination: { [key: string]: string; } = {
    "Transcriptome Profiling": "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\transcriptome profiling",
    "Biospecimen": "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\biospecimen",
    "Clinical":"C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\clinical",
    "Copy Number Variation": "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\copy number variation",
    "DNA Methylation": "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\dna methylation",
    "Simple Nucleotide Variation": "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\simple nucleotide variation"


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
        const {file_name, data_category } = record;

        
        
        const sourceFilePath: string = join(sourceFolderPath, file_name);

        await access(sourceFilePath, constants.R_OK);
        
        const destinationFolderPath: string = dataCategoryDestination[data_category];
        
        if(!destinationFolderPath) {
            throw new Error (`No destination folder for ${data_category}`);
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
        console.log(`Failed to process ${record.file_name}`)
    }

        

        
    }

})();