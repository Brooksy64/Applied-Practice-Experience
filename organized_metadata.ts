import { parse } from "path";

const{ promises: { readdir, csv-parse, readFile } } = require('fs');


const fileSourcePath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\dna methylation by platform\\HM450"
const samplesheetSourcePath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\gdc_sample_sheet.2022-06-08.csv"
const clinicalsheetSourcePath: string = "C:\\Users\\ericb\\Documents\\APE\\TARGET-AML\\data types\\clinical\\TARGET_AML_ClinicalData_AML1031_20211201.csv"

type TargetIDRecord = {
    file_name: string;
    case_id: string;

};

(async ()) => {
    //read filesystem in as a string
    const hm450FolderContents = await readdir(fileSourcePath, {withFileTypes: true});

    //read each file in as a string? array?
   

    //
}