import sys
from Bio import Entrez

# Set up Entrez with the email (required by NCBI)
Entrez.email = sys.argv[1]

# Get article IDs from command-line arguments
article_ids = sys.argv[2:]

# Function to fetch article abstracts
def fetch_abstracts(article_ids):
    abstracts = []
    for article_id in article_ids:
        handle = Entrez.esummary(db="pubmed", id=article_id)
        record = Entrez.read(handle)
        title = record[0].get("Title", "")
        source = record[0].get("Source", "")
        abstract = f"{title}: {source}"
        abstracts.append(abstract)
    return abstracts

# Fetch abstracts
abstracts = fetch_abstracts(article_ids)

# Output abstracts
print("\n".join(abstracts))
