# AWS And Nextflow Scaling

HEKBlueR can start as a local Shiny app and scale into a cloud screening pipeline.

## Local Mode

Use local mode for small runs.

```bash
nextflow run main.nf -profile local
```

## Docker Mode

Use Docker for reproducible dependencies.

```bash
docker build -t hekbluer:latest .
nextflow run main.nf -profile docker
```

## AWS Batch Mode

Use AWS Batch for many plates or scheduled processing.

Recommended services:

- S3 for raw uploads and result archives
- ECR for the HEKBlueR Docker image
- AWS Batch for analysis jobs
- Spot instances for compute tasks
- small on-demand compute for the Nextflow head job
- RDS Postgres for searchable run summaries
- CloudWatch for logs
- Secrets Manager for database credentials
- Cognito or SSO for secure app access

## Spot-Aware Design

Use Spot instances for batch analysis jobs because plate-level analysis is restartable and low memory.

Recommended settings:

- separate on-demand and Spot queues
- `aws.batch.maxSpotAttempts = 5`
- `process.errorStrategy = 'retry'`
- `process.maxRetries = 3`
- small CPU and memory requests per task
- S3 lifecycle rules for work directories

## Shiny Plus Batch Workflow

1. Biologist uploads files through Shiny.
2. Shiny writes raw data and metadata to S3.
3. Shiny writes run summary to Postgres.
4. Shiny starts a Nextflow or AWS Batch job.
5. AWS Batch runs the HEKBlueR pipeline.
6. Results are written to S3.
7. Key result tables are loaded into Postgres.
8. Shiny displays completed runs.

## Files To Edit Before AWS Use

Edit `nextflow.config`:

- AWS region
- S3 work directory
- AWS Batch queue
- ECR container image

The included AWS profiles are templates. They should not be used with placeholder bucket or container names.

