version 1.0

workflow pubmed_gpt3_literature_mining {
  input {
    String search_term
    Int max_results
    String openai_api_key
    String email
  }

  call fetch_pubmed_articles {
    input:
      search_term = search_term
      max_results = max_results
      email = email
  }

  call fetch_abstracts {
    input:
      article_ids = fetch_pubmed_articles.article_ids
  }

  call extract_disease_gene_associations {
    input:
      abstracts = fetch_abstracts.abstracts
      openai_api_key = openai_api_key
  }

  call save_results_to_csv {
    input:
      disease_gene_associations = extract_disease_gene_associations.results
  }

  output {
    File output_csv = save_results_to_csv.output_file
  }
}

task fetch_pubmed_articles {
  input {
    String search_term
    Int max_results
    String email
  }

  command <<<
    python3 fetch_pubmed_articles.py --search_term ${search_term} --max_results ${max_results} --email ${email}
  >>>

  output {
    Array[String] article_ids = read_lines(stdout())
  }
}

task fetch_abstracts {
  input {
    Array[String] article_ids
  }

  command <<<
    python3 fetch_abstracts.py --article_ids ${sep=" " article_ids}
  >>>

  output {
    Array[String] abstracts = read_lines(stdout())
  }
}

task extract_disease_gene_associations {
  input {
    Array[String] abstracts
    String openai_api_key
  }

  command <<<
    python3 extract_disease_gene_associations.py --abstracts ${sep=" " abstracts} --openai_api_key ${openai_api_key}
  >>>

  output {
    Array[String] results = read_lines(stdout())
  }
}

task save_results_to_csv {
  input {
    Array[String] disease_gene_associations
  }

  command <<<
    python3 save_results_to_csv.py --disease_gene_associations ${sep=" " disease_gene_associations}
  >>>

  output {
    File output_file = "disease_gene_associations.csv"
  }
}
