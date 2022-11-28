### Current metadata convention resides in schema/&lt;project&gt;_convention

* &lt;project&gt;_convention_schema.tsv - human-readable file used to generate the JSON document
* &lt;project&gt;_convention_schema.json - used by ingest pipeline for metadata validation
* scp_bq_inputs.json - file indicating non-metadata convention terms added to BigQuery for SCP use
* &lt;project&gt;_convention_schema.bq_schema.json - representation of BigQuery schema

# To update an existing metadata convention

* Create a new snapshot directory under schema/&lt;project&gt;_convention/snapshot using semantic versioning conventions  

* Make desired metadata convention updates to a copy of &lt;project&gt;_convention_schema.tsv file in the snapshot directory  

* copy scp_bq_inputs.json from previous snapshot directory, update with new SCP-internal terms if appropriate

* In the scripts `scp-ingest-pipeline/schema` directory, run

  ```
  python ../scripts/serialize_convention.py <project> <version>
  ```
  
* Copy the new convention JSON and TSV files to the * &lt;project&gt;_convention directory  
  
Notes:

* Tests in test_validate_metadata.py use current metadata convention (except for invalid metadata convention test)

* Specifically `test_bigquery_json_content` is expected to fail when the metadata convention is updated. The reference file, bq_test.json, must be updated (replace existing file with the generated addedfeed000000000000000.json file)

```
python validate_metadata.py --bq-json <path to metadata file>
```

* To create updated issues.json files to update reference files for tests, in the ingest/validation directory, run

```
python validate_metadata.py --issues-json <path to metadata file>
```

* To run validate_metadata.py against a different convention file:

```
python validate_metadata.py --convention <path to convention file> <path to metadata file>
```
