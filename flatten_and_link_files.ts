import { fstat, write } from "fs";

const { promises: { readFile, writeFile } } = require('fs');
const {join} = require('path');
const { parse } = require('csv/sync');

const bioFilePath: string =         "D:\\Documents 2\\APE\\TARGET-AML\\biospecimen.cart.2022-06-08.json";
const clinicalFilePath: string =    "D:\\Documents 2\\APE\\TARGET-AML\\clinical.cart.2022-06-08.json";
const sampleFilePath: string =      "D:\\Documents 2\\APE\\TARGET-AML\\gdc_sample_sheet.2022-06-08.csv";

const fileDestinationPath: string = "D:\\Documents 2\\APE\\TARGET-AML\\joined_record.json";
const fileDestinationPathMI: string = "D:\\Documents 2\\APE\\TARGET-AML\\joined_record_masked_intensities.json";
type bioRecord = {
    case_id: string;
    samples: ({
        sample_type_id: string;
        sample_id: string;
        freezing_method: string;
        submitter_id: string;
        tumor_code_id: string;
        sample_type: string;
        tumor_code: string;
        oct_embedded: string;
        updated_datetime: string;
        days_to_sample_procurement: string;
        days_to_collection: string;
        state: string;
        current_weight: string;
        is_ffpe: string;
        initial_weight: string;
        tissue_type: string;
    })[];
    submitter_id: string;
    sample_type: string;

};

type clinicalRecord = {
    case_id: string;
    diagnoses: ({
        primary_diagnosis: string;
    })[];
    demographic: {
        cause_of_death: string;
        race: string;
        gender: string;
        ethnicity: string;
        vital_status: string;
        age_at_index: string;
        submitter_id: string;
        days_to_birth: string;
        created_datetime: string;
        year_of_birth: string;
        premature_at_birth: string;
        weeks_gestation_at_birth: string;
        demographic_id: string;
        updated_datetime: string;
        days_to_death: string;
        state: string;
        year_of_death: string;
    };
};

type sampleRecord = {
    file_name: string;
    sample_id: string;
    data_category: string;
    data_type: string;
    project_id: string;
};


(async() => {
    //read json in as a string
    const bioFileContents = await readFile(bioFilePath, 'utf8');
    const clinicalFileContents = await readFile(clinicalFilePath, 'utf8');
    const sampleFileContents = await readFile(sampleFilePath, 'utf8');

    //parse the string into a JSON object
    try {
        const bio: bioRecord[] = JSON.parse(bioFileContents);
        const clinical: clinicalRecord[] = JSON.parse(clinicalFileContents);
        const sample: string[][] = parse(sampleFileContents);
        //console.log(sample)
        
        const bioSamplesMap = bio
             .flatMap(bioSamples => bioSamples.samples.map(s => ({ sample_type: s.sample_type, submitter_id: s.submitter_id, case_id: bioSamples.case_id})));
        //console.log(bioSamplesMap);

        const clinicalDemographicsMap = clinical
        .flatMap(clinicalSamples => clinicalSamples.diagnoses.map(d => ({primary_diagnosis: d.primary_diagnosis, race: clinicalSamples.demographic?.race, gender: clinicalSamples.demographic?.gender, ethnicity: clinicalSamples.demographic?.ethnicity, vital_status: clinicalSamples.demographic?.vital_status, case_id: clinicalSamples.case_id })));
        //console.log(clinicalDeomgraphicsMap)

        const serializedSamples: sampleRecord[] = sample.map( row => ({file_name: row[1], sample_id: row[6], data_category: row[2], data_type: row[3], project_id: row[4]}));
        //console.log(serializedSamples)

        const joinedRecord: any[] = serializedSamples.map(aVal => {
            const link1 = bioSamplesMap.find(bVal => bVal.submitter_id == aVal.sample_id.split(",")[0]);
            const link2 = clinicalDemographicsMap.find(cVal => cVal.case_id == link1?.case_id);

            return {
                ...aVal,
                ...link1,
                ...link2
            };
        });
        //console.log(joinedRecord);

        await writeFile(fileDestinationPath, JSON.stringify(joinedRecord));

        console.log('complete');

      const joinedRecordMaskedIntensities: any[] = joinedRecord.filter(m => m.data_type == "Masked Intensities" && !m.file_name.endsWith("_Red.idat") ) .map(e => {
e.file_name = e.file_name.replace("_Grn.idat","");

        return e;
      });
    
      
      await writeFile(fileDestinationPathMI, JSON.stringify(joinedRecordMaskedIntensities));

      console.log('complete');
        


    }
    catch(error) {
        console.log('error parsing JSON', error);
    }

    
     }
    

)();

