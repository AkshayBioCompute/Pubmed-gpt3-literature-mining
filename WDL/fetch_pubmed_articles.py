import sys
from Bio import Entrez

# Set up Entrez with the email (required by NCBI)
Entrez.email = sys.argv[3]

# Get search term and max results from command-line arguments
search_term = sys.argv[1]
max_results = int(sys.argv[2])

# Function to fetch PubMed article IDs
def fetch_pubmed_articles(search_term, max_results):
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=max_results)
    record = Entrez.read(handle)
    article_ids = record["IdList"]
    return article_ids

# Fetch articles
article_ids = fetch_pubmed_articles(search_term, max_results)

# Output article IDs
print("\n".join(article_ids))
