from .identification_pipelines.identification_pipeline_blastn import identification_pipeline_blastn

def run_identification(input_name: str, db: str = None):
    identification_pipeline_blastn(input_name, db)