library(shiny)
library(shinyjs)
library(DT)
library(stringr)
library(dplyr)
library(Biostrings)
library(bslib)

# ———————————————— PYTHON VENV SETUP ————————————————
# Creates and activates a Python virtualenv with numpy, biopython, tensorflow
venv_name   <- "r-reticulate"
venv_path   <- file.path("~", ".virtualenvs", venv_name)
system_python <- "/usr/bin/python3"

if (!dir.exists(venv_path)) {
  message("Creating Python venv at ", venv_path)
  reticulate::virtualenv_create(envname = venv_name, python = system_python)
  reticulate::virtualenv_install(
    envname  = venv_name,
    packages = c("pip", "wheel", "setuptools", "numpy", "biopython", "tensorflow"),
    ignore_installed = FALSE
  )
}
venv_python <- file.path(venv_path, "bin", "python")
Sys.setenv(RETICULATE_PYTHON = venv_python)
message("RETICULATE_PYTHON set to ", Sys.getenv("RETICULATE_PYTHON"))

library(reticulate)
use_python(venv_python, required = TRUE)
py_available(initialize = TRUE)


#------------------------------------------------------------------------------#
# 0. parse_mutation (verbatim, with fixed "ins" branch)
#------------------------------------------------------------------------------#
parse_mutation <- function(mutation, wildtype_seq) {
  info <- list(
    Locus            = NA_integer_,
    Mutation_Type    = NA_character_,
    New_Base         = NA_character_,
    Original_Base    = NA_character_,
    Deleted_Bases    = NA_character_,
    Duplicated_Bases = NA_character_,
    Inserted_Bases   = NA_character_,
    Mutated_Sequence = NA_character_
  )
  if (str_starts(mutation, "c.")) mutation <- substring(mutation, 3)
  if (str_detect(mutation, ">")) {
    parts <- str_split(mutation, ">")[[1]]
    if (length(parts) == 2) {
      left <- parts[1]; nb <- parts[2]
      num <- str_sub(left, 1, -2)
      if (str_detect(num, "^[0-9]+$") && str_detect(nb, "^[ACGT]$")) {
        loc <- as.integer(num)
        if (loc >= 1 && loc <= str_length(wildtype_seq)) {
          orig <- str_sub(wildtype_seq, loc, loc)
          mut  <- paste0(
            str_sub(wildtype_seq, 1, loc - 1),
            nb,
            str_sub(wildtype_seq, loc + 1)
          )
          info <- modifyList(info, list(
            Locus            = loc,
            Mutation_Type    = "Point",
            New_Base         = nb,
            Original_Base    = orig,
            Mutated_Sequence = mut
          ))
        }
      }
    }
  } else if (str_detect(mutation, "del") && !str_detect(mutation, "delins")) {
    parts <- str_split(mutation, "del")[[1]]
    rp    <- parts[1]; extra <- if (length(parts) > 1) parts[2] else ""
    if (str_detect(rp, "_")) {
      pos <- as.integer(str_split(rp, "_")[[1]])
      if (length(pos)==2 && !any(is.na(pos))) {
        db <- str_sub(wildtype_seq, pos[1], pos[2])
        mut <- paste0(
          str_sub(wildtype_seq, 1, pos[1]-1),
          str_sub(wildtype_seq, pos[2]+1)
        )
        info <- modifyList(info, list(
          Locus            = pos[1],
          Mutation_Type    = "Deletion",
          Deleted_Bases    = db,
          Mutated_Sequence = mut
        ))
      }
    } else {
      p <- as.integer(rp)
      if (!is.na(p)) {
        db <- if (extra != "" && str_detect(extra, "^[ACGT]$")) extra else str_sub(wildtype_seq, p, p)
        mut <- paste0(
          str_sub(wildtype_seq, 1, p-1),
          str_sub(wildtype_seq, p+1)
        )
        info <- modifyList(info, list(
          Locus            = p,
          Mutation_Type    = "Deletion",
          Deleted_Bases    = db,
          Mutated_Sequence = mut
        ))
      }
    }
  } else if (str_detect(mutation, "dup")) {
    parts <- str_split(mutation, "dup")[[1]]
    left  <- parts[1]
    if (str_detect(left, "_")) {
      pos <- as.integer(str_split(left, "_")[[1]])
      if (length(pos)==2 && !any(is.na(pos))) {
        db <- str_sub(wildtype_seq, pos[1], pos[2])
        mut <- paste0(
          str_sub(wildtype_seq, 1, pos[2]),
          db,
          str_sub(wildtype_seq, pos[2]+1)
        )
        info <- modifyList(info, list(
          Locus            = pos[1],
          Mutation_Type    = "Duplication",
          Duplicated_Bases = db,
          Mutated_Sequence = mut
        ))
      }
    } else {
      p <- as.integer(left)
      if (!is.na(p)) {
        db <- str_sub(wildtype_seq, p, p)
        mut <- paste0(
          str_sub(wildtype_seq, 1, p),
          db,
          str_sub(wildtype_seq, p+1)
        )
        info <- modifyList(info, list(
          Locus            = p,
          Mutation_Type    = "Duplication",
          Duplicated_Bases = db,
          Mutated_Sequence = mut
        ))
      }
    }
  } else if (str_detect(mutation, "delins")) {
    parts <- str_split(mutation, "delins")[[1]]
    rp    <- parts[1]; ins <- parts[2]
    pos   <- as.integer(str_split(rp, "_")[[1]])
    if (length(pos)==1) pos <- c(pos, pos)
    if (length(pos)==2 && !any(is.na(pos))) {
      db  <- str_sub(wildtype_seq, pos[1], pos[2])
      mut <- paste0(
        str_sub(wildtype_seq, 1, pos[1]-1),
        ins,
        str_sub(wildtype_seq, pos[2]+1)
      )
      info <- modifyList(info, list(
        Locus            = pos[1],
        Mutation_Type    = "Indel",
        Deleted_Bases    = db,
        Inserted_Bases   = ins,
        Mutated_Sequence = mut
      ))
    }
  } else if (str_detect(mutation, "ins") && !str_detect(mutation, "delins")) {
    parts <- str_split(mutation, "ins")[[1]]
    if (length(parts) == 2) {
      left <- parts[1]; ins <- parts[2]
      pos_vals <- as.integer(str_split(left, "_")[[1]])
      if (length(pos_vals) == 1 && !is.na(pos_vals)) {
        start <- pos_vals
        end   <- pos_vals
      } else if (length(pos_vals) == 2 && !any(is.na(pos_vals))) {
        start <- pos_vals[1]
        end   <- pos_vals[2]
      } else {
        return(info)
      }
      # Insert between start and end
      mut <- paste0(
        str_sub(wildtype_seq, 1, start),
        ins,
        str_sub(wildtype_seq, end + 1)
      )
      info <- modifyList(info, list(
        Locus            = start,
        Mutation_Type    = "Insertion",
        Inserted_Bases   = ins,
        Mutated_Sequence = mut
      ))
    }
  }
  info
}

#------------------------------------------------------------------------------#
# 1. translate_prot
#------------------------------------------------------------------------------#
translate_prot <- function(dna_seq) {
  if (is.na(dna_seq) || nchar(dna_seq) < 3) return(NA_character_)
  prot_full <- suppressWarnings(
    as.character(translate(DNAString(dna_seq), if.fuzzy.codon = "X"))
  )
  if (grepl("\\*", prot_full)) {
    cut <- regexpr("\\*", prot_full)[1] - 1
    substr(prot_full, 1, cut)
  } else prot_full
}

#------------------------------------------------------------------------------#
# 2. Refinement helpers
#------------------------------------------------------------------------------#
refine_point <- function(locus, seq, wt_seq, codon_tbl) {
  idx <- floor((locus - 1) / 3) * 3 + 1
  if (is.na(locus) || (idx + 2) > nchar(wt_seq) || (idx + 2) > nchar(seq))
    return("Point")
  wt_cod  <- toupper(substr(wt_seq, idx, idx + 2))
  mut_cod <- toupper(substr(seq,    idx, idx + 2))
  if (!wt_cod %in% names(codon_tbl) || !mut_cod %in% names(codon_tbl))
    return("Point")
  wt_aa  <- codon_tbl[[wt_cod]]
  mut_aa <- codon_tbl[[mut_cod]]
  if (wt_aa == mut_aa)      "Silent"
  else if (mut_aa == "*")   "Nonsense"
  else                       "Missense"
}

refine_indel <- function(mt, del, dup, ins) {
  net <- switch(mt,
                Deletion    = -nchar(del),
                Duplication =  nchar(dup),
                Insertion   =  nchar(ins),
                Indel       =  nchar(ins) - nchar(del),
                0)
  if (net %% 3 == 0) "In-Frame indel" else "Frameshifting indel"
}

#------------------------------------------------------------------------------#
# UI
#------------------------------------------------------------------------------#
ui <- fluidPage(
  useShinyjs(),
  theme = bs_theme(bootswatch = "darkly"),
  tags$head(
    tags$style(HTML("
      .shiny-progress .progress { height: 40px; margin-top: 10px; }
      .shiny-progress .progress-bar {
        background-image: linear-gradient(45deg, rgba(255,255,255,0.15) 25%, transparent 25%, transparent 50%, rgba(255,255,255,0.15) 50%, rgba(255,255,255,0.15) 75%, transparent 75%, transparent);
        background-size: 2rem 2rem;
        animation: progress-bar-stripes 1s linear infinite;
      }
      @keyframes progress-bar-stripes { from { background-position: 2rem 0; } to { background-position: 0 0;} }
    "))
  ),
  
  fluidRow(
    column(12, align = "center",
           tags$br(),
           tags$h1("MAP 🗺", style = "font-size:3em; margin-bottom:0; font-weight:bold;"),
           tags$h4("Molecular feature Association Pipeline️", style = "margin-top:0; color:#ccc; font-weight:bold;"),
           tags$div(
             "⚠️ Please cite: Schubert et al. (2025) [CITATION]",
             style = "font-size:1.1em; color:#fff; margin-top:5px; margin-bottom:5px;"
           ),
           tags$div(
             class = "alert alert-success alert-dismissible fade show",
             role  = "alert",
             style = "max-width:700px; margin:20px auto 0; padding:1rem;",
             tags$button(
               type             = "button",
               class            = "btn-close",
               `data-bs-dismiss` = "alert",
               `aria-label`     = "Close"
             ),
             HTML(paste0(
               "About 8.9% of the human genome are single exon genes",
               "<sup><a href='https://doi.org/10.1093/database/baw095' target='_blank' ",
               "title='Jorquera, R., Ortiz, R., Ossandon, F., Cárdenas, J. P., Sepúlveda, R., González, C., & Holmes, D. S. (2016). SinEx DB: a database for single exon coding sequences in mammalian genomes. Database: The Journal of Biological Databases and Curation, baw095.'>",
               "[1]",
               "</a></sup>",
               ". For those genes, we expect mRNA transcripts with premature stop codons not to be subject to nonsense-mediated decay but instead result in the presence of an altered, often truncated, protein",
               "<sup><a href='https://doi.org/10.1038/nsmb.1550' target='_blank' ",
               "title='Brogna, S., & Wen, J. (2009). Nonsense-mediated mRNA decay (NMD) mechanisms. Nature Structural & Molecular Biology, 16(2), 107–113.'>",
               "[2]",
               "</a>, ",
               "<a href='https://doi.org/10.1111/dmcn.16018' target='_blank' ",
               "title='Schubert, T., & Schaaf, C. P. (2025). MAGEL2 (patho-)physiology and Schaaf-Yang syndrome. Developmental Medicine and Child Neurology, 67(1), 35–48.'>",
               "[3]",
               "</a></sup>",
               ". Thus, information about the resulting protein product may be helpful for investigating genotype-phenotype relationships. MAP helps scientists generate this information."
             ))
           )
    )
  ),
  tags$br(),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$h4("Upload", style = "margin-top:0;"),
      
      # Dismissible info alert
      tags$div(
        class = "alert alert-info alert-dismissible fade show",
        role = "alert",
        tags$button(
          type = "button",
          class = "btn-close",
          `data-bs-dismiss` = "alert",
          `aria-label` = "Close"
        ),
        HTML(
          "<strong>Required files:<br><br></strong>
       <ul style='margin-bottom:0;'>
         <li>Variant CSV (with a <code>DNA_variant</code> column that contains variants in HGVS‐style cDNA notation, e.g. <code>c.123delA</code>)</li>
         <li>Reference FASTA (coding sequence)</li>
       </ul>
          <br><strong>You can try MAP with mock data:</strong>"
        ),
        # Example data buttons inside the same alert
        tags$div(
          style = "margin-top:10px;",
          fluidRow(
            column(6, actionButton("load_example_data", "Load example data", class = "btn-primary btn-block")),
            column(6, actionButton("clear_example_data", "Clear example data", class = "btn-secondary btn-block"))
          )
        )
      ),
      
      uiOutput("csv_ui"),
      uiOutput("fasta_ui"),
      checkboxInput("use_titer", "Include non-canonical TIS (TITER) analysis", FALSE),
      conditionalPanel(
        condition = "input.use_titer == true",
        uiOutput("fasta_flank_ui"),
        tags$br(),
        tags$br()
      ),
      actionButton("clear_uploads", "Clear Uploads", class = "btn-secondary"),
      tags$hr(),
      
      h4("Domains"),
      fluidRow(column(6, textInput("dom_name","Name")), column(6, numericInput("dom_start","Start (AA)",1,min=1))),
      fluidRow(column(6, numericInput("dom_end","End (AA)",1,min=1)), column(6, checkboxInput("dom_to_end","to end",FALSE))),
      fluidRow(column(6, actionButton("add_dom","Add",class="btn-warning btn-block")), column(6, actionButton("clr_dom","Clear",class="btn-secondary"))),
      tableOutput("dom_tbl") %>% tagAppendAttributes(style="font-size:0.9em; margin-top:10px; margin-bottom:20px;"),
      tags$hr(),
      
      h4("Motifs"),
      fluidRow(column(6, textInput("mot_name","Name")), column(6, textInput("mot_pattern","Pattern (AA seq.)"))),
      fluidRow(column(6, actionButton("add_mot","Add",class="btn-warning btn-block")), column(6, actionButton("clr_mot","Clear",class="btn-secondary"))),
      tableOutput("mot_tbl") %>% tagAppendAttributes(style="font-size:0.9em; margin-top:10px; margin-bottom:20px;"),
      tags$hr(),
      actionButton("run", "Run Analysis", icon = icon("play"), class = "btn-success btn-block", disabled = TRUE),
      actionButton("cancel", "Cancel Analysis", class = "btn-danger btn-block", disabled = TRUE)
    ),
    
    mainPanel(
      width = 9,
      DTOutput("res_tbl"),
      br(),
      uiOutput("download_btn")
    )
  ),
  tags$div(style="height:40px"),
  tags$div(
    style = "text-align:center; color:#ccc; margin-bottom:10px;",
    HTML(paste0(
      "This application uses the following R packages: shiny, shinyjs, DT, stringr, dplyr, Biostrings, bslib. ",
      "To predict non-canonical translation initiation sites, we deploy an adapted version of TITER ",
      "<sup><a href='https://doi.org/10.1093/bioinformatics/btx247' target='_blank' title='Zhang, S., Hu, H., Jiang, T., Zhang, L., & Zeng, J. (2017). TITER: predicting translation initiation sites by deep learning. Bioinformatics, 33(14), i234–i242.'>[4]</a>, </sup>",
      "<sup><a href='https://github.com/zhangsaithu/titer' target='_blank' title='TITER GitHub repository'>[5]</a></sup>",
      "."
    ))
  ),
  tags$div(style = "height:25px;"),
  tags$footer(
    style = "
    background: #fff;
    padding: 20px;
    display: flex;
    align-items: center;
    justify-content: space-between;
  ",
    # Left logo
    tags$div(
      style = "flex: 0 0 auto;",
      tags$img(
        src    = "https://lh3.googleusercontent.com/d/1-h0eXc4OnKe8dnS4xLrx8VvyAuN4eAX6",
        height = "80px"
      )
    ),
    
    # Center copyright + license
    tags$div(
      style = "
      flex: 1 1 auto;
      text-align: center;
      font-size: 0.7em;
      color: #555;
    ",
      # first paragraph
      tags$p(
        HTML('&copy; 2025 Tim Schubert.'),
        style = "margin: 0;"
      ),
      # second paragraph
      tags$p(
        tags$a(
          href   = "https://www.apache.org/licenses/LICENSE-2.0.html",
          target = "_blank",
          "Licensed under Apache License 2.0"
        ),
        style = "margin: 0;"
      )
    ),
    
    # Right logo
    tags$div(
      style = "flex: 0 0 auto;",
      tags$img(
        src    = "https://lh3.googleusercontent.com/d/1oXjxPtGKgD5TXdV8aHOFjjdhqr1e-OX5",
        height = "80px"
      )
    )
  )
)

#------------------------------------------------------------------------------#
# Server
#------------------------------------------------------------------------------#
server <- function(input, output, session) {
  # Reactive values to hold example data
  example_csv <- reactiveVal(NULL)
  example_fasta_seq <- reactiveVal(NULL)
  example_fasta_flank_seq <- reactiveVal(NULL)
  example_domains <- data.frame(Name = c("DomainA", "DomainB"), Start = c(5, 50), End = c(20, 100), stringsAsFactors = FALSE)
  example_motifs <- data.frame(Name = c("MotifX", "MotifY"), Pattern = c("ACDE", "WXYZ"), stringsAsFactors = FALSE)
  domains <- reactiveVal(data.frame(Name=character(),Start=integer(),End=numeric(),stringsAsFactors=FALSE))
  motifs  <- reactiveVal(data.frame(Name=character(),Pattern=character(),stringsAsFactors=FALSE))
  wt_aa_len <- reactiveVal(NULL)
  
  # Reactive to hold latest results
  result_data <- reactiveVal(NULL)
  cancel_requested <- reactiveVal(FALSE)
  observeEvent(input$cancel, {
    cancel_requested(TRUE)
    result_data(NULL)
    shinyjs::hide("res_tbl")
  })
  
  # Enable run button when required inputs are present (and flank FASTA if TITER is on)
  observe({
    has_csv <- !is.null(input$csv) || !is.null(example_csv())
    has_fasta <- !is.null(input$fasta) || !is.null(example_fasta_seq())
    has_flank <- if (isTRUE(input$use_titer)) {
      !is.null(input$fasta_flank) || !is.null(example_fasta_flank_seq())
    } else TRUE
    shinyjs::toggleState("run", has_csv && has_fasta && has_flank)
  })
  
  # Clear uploads
  observeEvent(input$clear_uploads, {
    shinyjs::reset("csv")
    shinyjs::reset("fasta")
    shinyjs::reset("fasta_flank")
    example_csv(NULL)
    example_fasta_seq(NULL)
    example_fasta_flank_seq(NULL)
    result_data(NULL)
    shinyjs::disable("run")
    shinyjs::disable("cancel")
  })
  
  # Load all example data (variants, reference, domains, motifs)
  observeEvent(input$load_example_data, {
    # Load example variants
    example_csv(read.csv("www/example_variants.csv", stringsAsFactors = FALSE))
    # Load example reference
    lines <- readLines("www/example_reference.fasta", warn = FALSE)
    example_fasta_seq(toupper(paste(lines[!startsWith(lines, ">")], collapse = "")))
    # Load example reference with flanks for TITER
    lines_flank <- readLines("www/example_reference_flank.fasta", warn = FALSE)
    example_fasta_flank_seq(toupper(paste(lines_flank[!startsWith(lines_flank, ">")], collapse = "")))
    # Load example domains & motifs
    domains(example_domains)
    motifs(example_motifs)
    showNotification("Example data loaded.", type = "message")
  })
  
  # Domains: bounds & duplicate check
  observeEvent(input$add_dom, {
    req(input$dom_name, input$dom_start)
    max_aa <- wt_aa_len()
    end_val <- if (input$dom_to_end) max_aa else input$dom_end
    if (!is.null(max_aa) && end_val > max_aa) {
      showNotification(paste0("Domain end (", end_val, ") exceeds protein length (", max_aa, ")."), type = "error")
      return()
    }
    df <- domains()
    if (any(df$Name == input$dom_name & df$Start == input$dom_start & df$End == end_val)) {
      showNotification("Exact duplicate domain exists.", type="error")
      return()
    }
    domains(rbind(df, data.frame(Name=input$dom_name, Start=input$dom_start, End=end_val, stringsAsFactors=FALSE)))
  })
  observeEvent(input$clr_dom, domains(data.frame(Name=character(),Start=integer(),End=numeric(),stringsAsFactors=FALSE)))
  output$dom_tbl <- renderTable(domains(), striped=TRUE, hover=TRUE)
  
  # Motifs: IUPAC AA & duplicate check
  observeEvent(input$add_mot, {
    req(input$mot_name, input$mot_pattern)
    # Only allow IUPAC AA letters: A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,B,Z,X
    if (!str_detect(input$mot_pattern, "^[ACDEFGHIKLMNPQRSTVWYBZX]+$")) {
      showNotification("Pattern must use only IUPAC amino acid letters A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,B,Z,X", type = "error")
      return()
    }
    df <- motifs()
    if (any(df$Name == input$mot_name & df$Pattern == input$mot_pattern)) {
      showNotification("Exact duplicate motif exists.", type="error")
      return()
    }
    motifs(rbind(df, data.frame(Name=input$mot_name, Pattern=input$mot_pattern, stringsAsFactors=FALSE)))
  })
  observeEvent(input$clr_mot, motifs(data.frame(Name=character(),Pattern=character(),stringsAsFactors=FALSE)))
  output$mot_tbl <- renderTable(motifs(), striped=TRUE, hover=TRUE)
  
  # Clear example data
  observeEvent(input$clear_example_data, {
    example_csv(NULL)
    example_fasta_seq(NULL)
    example_fasta_flank_seq(NULL)
    domains(data.frame(Name=character(),Start=integer(),End=numeric(),stringsAsFactors=FALSE))
    motifs(data.frame(Name=character(),Pattern=character(),stringsAsFactors=FALSE))
    result_data(NULL)
    showNotification("Example data cleared.", type = "message")
  })
  
  # (alert dismiss observer removed)
  
  # Main analysis: run and store results
  observeEvent(input$run, {
    cancel_requested(FALSE)
    shinyjs::disable("run")
    shinyjs::enable("cancel")
    shinyjs::hide("res_tbl")
    # Require either uploaded or example data
    req(!is.null(input$csv) || !is.null(example_csv()))
    
    # CSV checks
    if (!is.null(input$csv)) {
      vd <- read.csv(input$csv$datapath, stringsAsFactors=FALSE)
    } else {
      vd <- example_csv()
    }
    if (nrow(vd) == 0) {
      showModal(modalDialog(
        title = "Empty CSV",
        "Your variant list file has zero rows. Please upload a non‐empty CSV.",
        easyClose = TRUE, footer = modalButton("OK")
      ))
      result_data(NULL)
      return(NULL)
    }
    if (!"DNA_variant" %in% colnames(vd)) {
      showModal(modalDialog(
        title = "Missing column",
        "Your CSV must contain a column named `DNA_variant`. Please fix and re‐upload.",
        easyClose = TRUE, footer = modalButton("OK")
      ))
      result_data(NULL)
      return(NULL)
    }
    if (isTRUE(input$use_titer)) {
      if (!"patient_ID" %in% colnames(vd)) {
        showModal(modalDialog(
          title = "Missing column",
          "TITER analysis requires a `patient_ID` column in your CSV.",
          easyClose = TRUE, footer = modalButton("OK")
        ))
        result_data(NULL); return(NULL)
      }
      if (is.null(input$fasta_flank) && is.null(example_fasta_flank_seq())) {
        showModal(modalDialog(
          title = "Missing FASTA",
          "Please upload the FASTA with 100 bp flanks for TITER analysis.",
          easyClose = TRUE, footer = modalButton("OK")
        ))
        result_data(NULL); return(NULL)
      }
    }
    
    # FASTA reading
    if (!is.null(input$fasta)) {
      lines <- readLines(input$fasta$datapath, warn = FALSE)
    } else {
      lines <- c(paste0(">example"), example_fasta_seq())
    }
    if (!any(startsWith(lines, ">"))) {
      showModal(modalDialog(
        title = "Invalid FASTA",
        "Your FASTA must have at least one header line beginning with `>`.",
        easyClose = TRUE, footer = modalButton("OK")
      ))
      result_data(NULL)
      return(NULL)
    }
    seq_nt <- toupper(paste(lines[!startsWith(lines, ">")], collapse = ""))
    if (nchar(seq_nt) %% 3 != 0) {
      showModal(modalDialog(
        title = "Sequence not multiple of 3",
        "Reference sequence length is not divisible by 3—please check your FASTA.",
        easyClose = TRUE, footer = modalButton("OK")
      ))
      result_data(NULL)
      return(NULL)
    }
    aa_len <- nchar(translate(DNAString(seq_nt), if.fuzzy.codon="X"))
    wt_aa_len(aa_len)
    
    withProgress(message="Running analysis...", value=0, {
      incProgress(0.1, detail="Parsing mutations")
      if (cancel_requested()) {
        showNotification("Analysis cancelled.", type = "warning")
        shinyjs::disable("cancel")
        result_data(NULL)
        return(NULL)
      }
      parsed <- lapply(vd$DNA_variant, parse_mutation, wildtype_seq = seq_nt)
      df <- vd %>% mutate(
        Locus = sapply(parsed, `[[`,"Locus"),
        Mutation_Type = sapply(parsed, `[[`,"Mutation_Type"),
        Deleted_Bases = sapply(parsed, `[[`,"Deleted_Bases"),
        Duplicated_Bases = sapply(parsed, `[[`,"Duplicated_Bases"),
        Inserted_Bases = sapply(parsed, `[[`,"Inserted_Bases"),
        Mutated_Sequence = sapply(parsed, `[[`,"Mutated_Sequence")
      )
      
      incProgress(0.2, detail="Translating proteins")
      if (cancel_requested()) {
        showNotification("Analysis cancelled.", type = "warning")
        shinyjs::disable("cancel")
        result_data(NULL)
        return(NULL)
      }
      df <- df %>% mutate(
        Mutated_Protein   = sapply(Mutated_Sequence, translate_prot),
        Protein_Length_aa = nchar(Mutated_Protein)
      )
      
      incProgress(0.2, detail="Refining mutation types")
      if (cancel_requested()) {
        showNotification("Analysis cancelled.", type = "warning")
        shinyjs::disable("cancel")
        result_data(NULL)
        return(NULL)
      }
      df <- df %>% rowwise() %>% mutate(Mutation_Type = case_when(
        Mutation_Type == "Point"         ~ refine_point(Locus, Mutated_Sequence, seq_nt, GENETIC_CODE),
        Mutation_Type %in% c("Deletion","Duplication","Insertion","Indel") ~ refine_indel(Mutation_Type, Deleted_Bases, Duplicated_Bases, Inserted_Bases),
        TRUE                             ~ Mutation_Type
      )) %>% ungroup()
      
      incProgress(0.2, detail="Annotating domains & motifs")
      if (cancel_requested()) {
        showNotification("Analysis cancelled.", type = "warning")
        shinyjs::disable("cancel")
        result_data(NULL)
        return(NULL)
      }
      doms <- domains(); mot <- motifs()
      df <- df %>% mutate(
        AA_Position = floor((Locus - 1) / 3) + 1,
        Domain_Location_Of_Variant = mapply(function(p){
          hits <- doms$Name[p >= doms$Start & p <= doms$End]
          if (!length(hits)) "" else paste(hits, collapse='; ')
        }, AA_Position),
        Lost_Functional_Domains = mapply(function(p, mt){
          hits <- doms$Name[p >= doms$Start & p <= doms$End]
          if (!length(hits)) return(NA_character_)
          mt <- tolower(mt)
          if (mt %in% c('nonsense','frameshifting indel')){
            i <- which(doms$Name %in% hits); lost <- doms$Name[i:length(doms$Name)]
          } else if (mt %in% c('missense','in-frame indel')) lost <- hits else return(NA_character_)
          paste(lost, collapse='; ')
        }, AA_Position, Mutation_Type)
      )
      
      if (nrow(mot) > 0) {
        for (i in seq_len(nrow(mot))) {
          nm <- mot$Name[i]; pat <- mot$Pattern[i]
          df[[nm]] <- ifelse(is.na(df$Mutated_Protein), "", ifelse(str_detect(df$Mutated_Protein, pat), "1", "0"))
        }
        df <- df %>% rowwise() %>% mutate(`Lost Motifs` = {
          vals <- c_across(all_of(mot$Name))
          miss <- mot$Name[vals == "0"]
          if (!length(miss)) "" else paste(miss, collapse='; ')
        }) %>% ungroup()
      } else {
        df$`Lost Motifs` <- ""
      }
      
      # --- Non-canonical TIS via TITER ---
      if (isTRUE(input$use_titer)) {
        incProgress(0.05, detail="Running non-canonical TIS analysis ⚠️ do not leave this page! TITER will take about 2 min for 10 patients)")
        if (cancel_requested()) {
          showNotification("Analysis cancelled.", type = "warning")
          shinyjs::disable("cancel")
          result_data(NULL)
          return(NULL)
        }
        # Copy titer folder to a temp directory so we can write files on shinyapps.io
        temp_titer_dir <- file.path(tempdir(), "titer")
        if (dir.exists(temp_titer_dir)) unlink(temp_titer_dir, recursive=TRUE)
        file.copy("titer", tempdir(), recursive=TRUE)
        # Ensure the data subdirectory exists before writing files
        dir.create(file.path(temp_titer_dir, "data"), recursive = TRUE, showWarnings = FALSE)
        # Copy uploaded or example flanked FASTA into data/ (preserve original filename)
        if (!is.null(input$fasta_flank)) {
          flank_basename <- basename(input$fasta_flank$name)
          file.copy(
            input$fasta_flank$datapath,
            file.path(temp_titer_dir, "data", flank_basename),
            overwrite = TRUE
          )
        } else {
          flank_basename <- "example_reference_flank.fasta"
          file.copy(
            "www/example_reference_flank.fasta",
            file.path(temp_titer_dir, "data", flank_basename),
            overwrite = TRUE
          )
        }
        # Write out only the columns needed by TITER
        write.csv(
          df %>% dplyr::select(patient_ID, Mutated_Sequence, DNA_variant),
          file.path(temp_titer_dir, "data", "variant_list_with_mutated_sequences.csv"),
          row.names = FALSE
        )
        # Run the Python TITER script
        script_path <- file.path(temp_titer_dir, "codes", "analyze_patients_for_variant_specific_additional_TIS.py")
        # Ensure Python environment is available
        if (!reticulate::py_available(initialize = FALSE)) {
          showNotification("Python not available; skipping TITER analysis.", type = "error")
          # Skip the TITER block by returning early or setting a flag
        } else {
          tryCatch({
            reticulate::use_virtualenv("titer-venv", required = TRUE)
          }, error = function(e) {
            message("Could not activate virtualenv; using existing Python environment. Reason: ", e$message)
          })
          reticulate::source_python(script_path)
          # Merge back the summary
          titer_summary <- read.csv(
            file.path(temp_titer_dir, "data", "summary_patients_most_likely_additional_TIS.csv"),
            stringsAsFactors = FALSE
          ) %>% dplyr::select(-DNA_variant)
          df <- left_join(df, titer_summary, by = "patient_ID")
          # Translate and add metrics for non-canonical TIS proteins
          df <- df %>% mutate(
            Protein_from_most_likely_non_canonical_TIS         = sapply(RNA_sequence_most_likely_non_canonical_TIS, translate_prot),
            Protein_from_most_likely_non_canonical_in_frame_TIS = sapply(RNA_sequence_most_likely_non_canonical_in_frame_TIS, translate_prot),
            Protein_Length_from_most_likely_non_canonical_TIS         = nchar(Protein_from_most_likely_non_canonical_TIS),
            Protein_Length_from_most_likely_non_canonical_in_frame_TIS = nchar(Protein_from_most_likely_non_canonical_in_frame_TIS)
          ) %>% {
            # Compute frameshifts and lost domains exactly as in the .Rmd
            df_inner <- .
            df_inner <- df_inner %>% mutate(
              Mutation_AA_Pos_Canonical = floor((Locus - 1) / 3) + 1,
              Main_TIS_bp    = as.numeric(Most_likely_non_canonical_TIS_CDS_Position),
              Main_TIS_Codon = floor((Main_TIS_bp - 1) / 3) + 1,
              Mutation_AA_Pos_Main_TIS = if_else(
                Mutation_Type == "Frameshifting indel" &
                  Mutation_AA_Pos_Canonical >= Main_TIS_Codon,
                Mutation_AA_Pos_Canonical - (Main_TIS_Codon - 1),
                NA_integer_
              ),
              Unaltered_Raw_Main  = Mutation_AA_Pos_Main_TIS - 1,
              Frameshift_Raw_Main = Protein_Length_from_most_likely_non_canonical_TIS - Unaltered_Raw_Main,
              Unaltered_Length_from_most_likely_non_canonical_TIS = if_else(
                Mutation_Type == "Frameshifting indel" &
                  !is.na(Unaltered_Raw_Main) &
                  Unaltered_Raw_Main >= 0 &
                  Unaltered_Raw_Main <= Protein_Length_from_most_likely_non_canonical_TIS,
                Unaltered_Raw_Main,
                NA_integer_
              ),
              Frameshift_Length_from_most_likely_non_canonical_TIS = if_else(
                Mutation_Type == "Frameshifting indel" &
                  !is.na(Frameshift_Raw_Main) &
                  Frameshift_Raw_Main >= 0,
                Frameshift_Raw_Main,
                NA_integer_
              )
            ) %>% mutate(
              InFrame_TIS_bp    = as.numeric(Most_likely_non_canonical_in_frame_TIS_CDS_Position),
              InFrame_TIS_Codon = floor((InFrame_TIS_bp - 1) / 3) + 1,
              Mutation_AA_Pos_Alternative = if_else(
                Mutation_Type == "Frameshifting indel" &
                  Mutation_AA_Pos_Canonical >= InFrame_TIS_Codon,
                Mutation_AA_Pos_Canonical - (InFrame_TIS_Codon - 1),
                NA_integer_
              ),
              Unaltered_Raw_InFrame  = Mutation_AA_Pos_Alternative - 1,
              Frameshift_Raw_InFrame = Protein_Length_from_most_likely_non_canonical_in_frame_TIS -
                (Unaltered_Raw_InFrame),
              Unaltered_Length_from_most_likely_non_canonical_in_frame_TIS = if_else(
                Mutation_Type == "Frameshifting indel" &
                  !is.na(Unaltered_Raw_InFrame) &
                  Unaltered_Raw_InFrame >= 0 &
                  Unaltered_Raw_InFrame <= Protein_Length_from_most_likely_non_canonical_in_frame_TIS,
                Unaltered_Raw_InFrame,
                NA_integer_
              ),
              Frameshift_Length_from_most_likely_non_canonical_in_frame_TIS = if_else(
                Mutation_Type == "Frameshifting indel" &
                  !is.na(Frameshift_Raw_InFrame) &
                  Frameshift_Raw_InFrame >= 0,
                Frameshift_Raw_InFrame,
                NA_integer_
              )
            )
            df_inner <- df_inner %>% dplyr::select(
              -Mutation_AA_Pos_Canonical, -Main_TIS_bp, -Main_TIS_Codon,
              -Mutation_AA_Pos_Main_TIS, -Unaltered_Raw_Main, -Frameshift_Raw_Main,
              -InFrame_TIS_bp, -InFrame_TIS_Codon, -Mutation_AA_Pos_Alternative,
              -Unaltered_Raw_InFrame, -Frameshift_Raw_InFrame
            )
            df_inner
          }
        }
      }
      
      incProgress(0.1, detail="Calculating frameshift metrics")
      if (cancel_requested()) {
        showNotification("Analysis cancelled.", type = "warning")
        shinyjs::disable("cancel")
        result_data(NULL)
        return(NULL)
      }
      df <- df %>% mutate(
        Unaltered_Length_aa  = if_else(Mutation_Type == 'Frameshifting indel', (floor((Locus-1)/3)+1)-1, NA_integer_),
        Frameshift_Length_aa = if_else(Mutation_Type == 'Frameshifting indel', Protein_Length_aa - ((floor((Locus-1)/3)+1)-1), NA_integer_)
      )
      
      colnames(df) <- gsub("_", " ", colnames(df))
      df <- df %>% select(-c(`Mutated Sequence`,`Mutated Protein`), everything(), `Mutated Sequence`,`Mutated Protein`)
      # At end of analysis, assign to result_data:
      result_data(df)
      shinyjs::disable("cancel")
      shinyjs::enable("run")
      shinyjs::show("res_tbl")
    })
  })
  
  output$res_tbl <- renderDT({
    req(result_data()); datatable(result_data(), rownames=FALSE, filter='top', options=list(pageLength=10, scrollX=TRUE))
  })
  
  output$dl <- downloadHandler(
    filename = function() paste0('MAP_results_', Sys.Date(), '.csv'),
    content  = function(f) write.csv(result_data(), f, row.names=FALSE)
  )
  
  output$download_btn <- renderUI({
    req(result_data())
    downloadButton("dl", "Download CSV", class = "btn-lg btn-success")
  })
  # UI outputs for file inputs, switching to "using example" if example data is loaded
  output$csv_ui <- renderUI({
    if (is.null(example_csv())) {
      fileInput(
        "csv",
        tagList(
          "Variant CSV ",
          tags$span(
            icon("info-circle", lib = "font-awesome"),
            title = "Upload a .csv file with a column named \"patient_ID\" (with individual IDs for each row) and another column \"DNA_variant\" containing variants in HGVS cDNA notation, e.g. \"c.123delA\".",
            style = "cursor: help;"
          )
        ),
        accept = c('.csv','text/csv')
      )
    } else {
      tags$div(
        style = "color:#ccc; font-size:0.9em; margin-bottom:10px;",
        "Using example Variant CSV"
      )
    }
  })
  
  output$fasta_ui <- renderUI({
    if (is.null(example_fasta_seq())) {
      fileInput(
        "fasta",
        tagList(
          "Reference FASTA (cds) ",
          tags$span(
            icon("info-circle", lib = "font-awesome"),
            title = paste(
              "Upload a FASTA file containing the coding sequence (CDS) of your reference gene.",
              "",
              "A typical approach can be:",
              "- find your gene sequence on ensembl.org",
              "- click download sequence",
              "- for included sequences: select all",
              "- select \"FASTA\" as the file format",
              "- set flanks to \"0\" in the dialogue",
              "- click \"Download\"",
              sep = "\n"
            ),
            style = "cursor: help;"
          )
        ),
        accept = c('.fa','.fasta','text/plain')
      )
    } else {
      tags$div(
        style = "color:#ccc; font-size:0.9em; margin-bottom:10px;",
        "Using example Reference FASTA (coding sequence)"
      )
    }
  })
  
  output$fasta_flank_ui <- renderUI({
    req(input$use_titer)
    if (is.null(example_fasta_flank_seq())) {
      tagList(
        fileInput(
          "fasta_flank",
          tagList(
            "Reference FASTA with 100 bp flanks ",
            tags$span(
              icon("info-circle", lib = "font-awesome"),
              title = paste(
                "Upload a FASTA file containing the coding sequence (CDS) of your reference gene with 100 bp flanks on either side.",
                "",
                "- find your gene sequence on ensembl.org",
                "- click download sequence",
                "- for included sequences: select all",
                "- select \"FASTA\" as the file format",
                "- set flanks to \"100\" in the dialogue",
                "- click \"Download\"",
                sep = "\n"
              ),
              style = "cursor: help;"
            )
          ),
          accept = c('.fa','.fasta','text/plain')
        ),
        tags$div(
          style = "font-size:0.8em; color:#ccc; margin-top:5px; margin-bottom:10px;",
          HTML(paste0(
            "Upload FASTA containing the coding sequence plus 100 bp upstream and downstream flanks for non-canonical TIS analysis. ⚠️ Please only run TITER when necessary, it will use extensive resources and take time. Do not leave this page during a run.",
            " (adapted from TITER ",
            "<sup><a href='https://doi.org/10.1093/bioinformatics/btx247' target='_blank'>[4]</a></sup>",
            "<sup><a href='https://github.com/zhangsaithu/titer' target='_blank'>[5]</a></sup>",
            ")"
          ))
        )
      )
    } else {
      tags$div(
        style = "color:#ccc; font-size:0.8em; margin-bottom:10px;",
        "Using example FASTA with 100 bp flanks for TITER. ⚠️ Please only run TITER when necessary, it will use extensive resources and take time. Do not leave this page during a run."
      )
    }
  })
}

shinyApp(ui, server)

