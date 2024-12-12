# Import necessary libraries
import openai
from Bio import Entrez
import csv

# Set your email for NCBI Entrez (required)
Entrez.email = "your_email@example.com"

# Set your OpenAI API key (replace with your actual key)
openai.api_key = 'your-openai-api-key'

# Function to fetch PubMed article IDs using a search term
def fetch_pubmed_articles(search_term, max_results=10):
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=max_results)
    record = Entrez.read(handle)
    article_ids = record["IdList"]
    return article_ids

# Function to fetch article abstracts by PubMed article IDs
def fetch_abstracts(article_ids):
    abstracts = []
    for article_id in article_ids:
        handle = Entrez.esummary(db="pubmed", id=article_id)
        record = Entrez.read(handle)
        abstract = record[0].get("Title", "") + ": " + record[0].get("Source", "")
        abstracts.append(abstract)
    return abstracts

# Function to use GPT-3 to extract disease-gene associations from abstracts
def extract_disease_gene_associations(abstracts):
    results = []
    for abstract in abstracts:
        response = openai.Completion.create(
            engine="text-davinci-003",  # Use GPT-3 engine
            prompt=f"Extract the disease-gene associations from the following scientific abstract: {abstract}",
            max_tokens=150,
            temperature=0.7
        )
        result = response.choices[0].text.strip()
        results.append(result)
    return results

# Function to post-process and save the results into a CSV file
def save_results_to_csv(disease_gene_associations, output_file="disease_gene_associations.csv"):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Abstract", "Disease-Gene Associations"])
        for abstract, associations in zip(disease_gene_associations[0], disease_gene_associations[1]):
            writer.writerow([abstract, associations])

# Main function to execute the workflow
def main():
    search_term = "Coronary Artery Disease"  # Example search term
    max_results = 10  # Number of articles to fetch
    
    # Step 1: Fetch PubMed article IDs
    article_ids = fetch_pubmed_articles(search_term, max_results)
    print(f"Fetched {len(article_ids)} article IDs.")
    
    # Step 2: Fetch article abstracts using the article IDs
    abstracts = fetch_abstracts(article_ids)
    print(f"Fetched {len(abstracts)} abstracts.")
    
    # Step 3: Extract disease-gene associations using GPT-3
    disease_gene_associations = extract_disease_gene_associations(abstracts)
    print(f"Extracted disease-gene associations from {len(disease_gene_associations)} abstracts.")
    
    # Step 4: Save the results to a CSV file
    save_results_to_csv(zip(abstracts, disease_gene_associations))
    print(f"Results saved to 'disease_gene_associations.csv'.")

# Run the main function to execute the workflow
if __name__ == "__main__":
    main()
