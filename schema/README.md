### Current metadata convention resides in schema/&lt;project&gt;_convention

* &lt;project&gt;_convention_schema.tsv - human-readable file used to generate the JSON document
* &lt;project&gt;_convention_schema.json - used by ingest pipeline for metadata validation

# To update an existing metadata convention

* Create a new snapshot directory under schema/&lt;project&gt;_convention/snapshot using semantic versioning conventions  

* Make desired metadata convention updates to a copy of &lt;project&gt;_convention_schema.tsv file in the snapshot directory  

* In the scripts `scp-ingest-pipeline` directory, run
  ```
  python serialize_convention.py <project> <version>
  ```
  
* Copy the new convention JSON and TSV files to the * &lt;project&gt;_convention directory  
  

* Copy new metadata schema to appropriate EXTERNAL_JSON_CONVENTION location 
```
gs://broad-singlecellportal-public/schema/<project>_convention/<version>/<project>_convention_schema.json
```  
  
* Update EXTERNAL_JSON_CONVENTION in ingest_pipeline.py

Notes:
* Tests in test_validate_metadata.py use current metadata convention (except for invalid metadata convention test)

* To create updated issues.json files to update reference files for tests, in the ingest/validation directory, run
```
python validate_metadata.py --issues-json <path to metadata file>
```

* To run validate_metadata.py against a different convention file:
```
python validate_metadata.py --convention <path to convention file> <path to metadata file>
```
