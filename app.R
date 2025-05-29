
library(shiny)
library(DT)
library(data.table)
library(VariantAnnotation)
library(ranger)
library(mlr3)
library(mlr3learners)
library(plotly)
library(umap)  # if not installed, auto-installed inside code
library(GenomicFeatures)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(enrichR)

options(shiny.maxRequestSize = 50000 * 1024^2)  # 500 MB
# SNP positions from final hits




ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),tags$head(
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Montserrat:wght@400;600;700&display=swap")
  ),
  
  
  titlePanel("JCAP GWAS Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("vcfFile", "Upload VCF File", accept = ".vcf"),
      fileInput("phenoFile", "Upload Phenotype Table (CSV)", accept = ".csv"),
      fileInput("covarFile", "Upload Covariates Table ", accept = ".csv"),
      numericInput("maf_thresh", "MAF Threshold", value = 0.01),
      
      numericInput("af_thresh", "ALT Allele Frequency Threshold", value = 0.01),      # NEW
      numericInput("hwe_thresh", "HWE P-value Threshold", value = 1e-6),              # NEW
      actionButton("runQC", "Run QC"),
      actionButton("runPCA", "Run PCA"),
      actionButton("runUMAP", "Run UMAP"),
      actionButton("runPower", "Run Power Analysis"),
      tabPanel("Power Analysis", DTOutput("powerTable")),
      
      actionButton("runPowerCurve", "Run Power Curve"),
      
      
      actionButton("runGWAS", "Run GWAS"),
      sliderInput("pval_thresh", "P-value Threshold",
                  min = 1e-10, max = 0.1, value = 0.05, step = 0.0001),
      
      checkboxInput("use_bonferroni", "Apply Bonferroni Correction", value = FALSE),
      
      actionButton("runQQ", "Run Q-Q Plot"),
      actionButton("runManhattan", "Run Manhattan Plot"),
      
      selectInput("chr_filter", "Chromosome", choices = c("All", as.character(1:22)), selected = "All"),
      numericInput("pos_min", "Min Position", value = 0),
      numericInput("pos_max", "Max Position", value = 1e6),
      actionButton("mapGenes", "Map SNPs to Nearest Genes"),
      actionButton("plotGenes", "Plot Gene Map"),
      actionButton("runEnrichment", "Run Enrichment"),
      actionButton("plotEnrichment", "Plot Enrichment Results"),
      selectInput(
        inputId = "enrich_db",
        label = "Choose enrichment database",
        choices = c("KEGG_2021_Human", "GO_Biological_Process_2021", "Reactome_2016"),
        selected = "Reactome_2016"),
      actionButton("runRF", "Train Random Forest"),
      actionButton("runROC", "Run ROC Plot"),
      
      
      downloadButton("download_power", "Download Power Table"),
      
      downloadButton("download_gwas", "Download GWAS Results"),
      downloadButton("download_covar", "Download Covariate P-Values"),
      downloadButton("download_final", "Download Final Hits"),
      downloadButton("download_annotated", "Download Annotated Hits"),
      downloadButton("download_enrichment", "Download Enrichment Table"),
      downloadButton("download_rf_preds", "Download RF Predictions"),
      downloadButton("download_rf_metrics", "Download RF Metrics"),
      downloadButton("download_rf_importance", "Download RF Importance"),
      
      
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(title = "README", value = "readme_tab", class = "readme-tab",
                 fluidRow(
                   column(
                     width = 12,
                     h2("📘 How to Use the JCAP GWAS App"),
                     tags$ol(
                       tags$li(strong("Upload Data:"), " Use the top three file inputs to upload your VCF file, phenotype table (CSV), and optional covariates table."),
                       tags$li(strong("Perform QC:"), " Set the MAF, AF, and HWE thresholds, then click ", code("Run QC"), " to filter your SNPs."),
                       tags$li(strong("Exploratory Visualization:"), " Click ", code("Run PCA"), " and ", code("Run UMAP"), " to view dimensionality-reduced sample clustering."),
                       tags$li(strong("Power Analysis:"), " Use ", code("Run Power Analysis"), " to estimate power based on a top SNP, and ", code("Run Power Curve"), " to visualize power vs. sample size."),
                       tags$li(strong("GWAS Execution:"), " Set your desired p-value threshold using the slider and (optionally) toggle Bonferroni correction, then click ", code("Run GWAS"), ". This will generate association p-values for all SNPs."),
                       tags$li(strong("Interpret Final Hits:"), " The ", strong("Final Phenotype Hits"), " table contains SNPs associated with the phenotype but not with any covariates."),
                       tags$li(strong("GWAS Visualizations:"), " Use ", code("Run Q-Q Plot"), " and ", code("Run Manhattan Plot"), " to view GWAS plots. You can filter the Manhattan plot by chromosome and base pair position."),
                       tags$li(strong("SNP-to-Gene Mapping:"), " Click ", code("Map SNPs to Nearest Genes"), " to associate significant SNPs with genes, and ", code("Plot Gene Map"), " to visualize them."),
                       tags$li(strong("Pathway Enrichment:"), " Click ", code("Run Enrichment"), " to perform pathway analysis, then choose a database and click ", code("Plot Enrichment Results"), " to visualize the top terms."),
                       tags$li(strong("Machine Learning (Random Forest):"), " Click ", code("Train Random Forest"), " to build a predictive model using final hits. This will output predicted classes, AUC, sensitivity/specificity, and SNP importance rankings. You can also generate an ROC curve with ", code("Run ROC Plot"), "."),
                       tags$li(strong("Download Results:"), " Use the download buttons on the sidebar and inside each tab to export all tables as CSVs.")
                     ),
                     tags$hr(),
                     p("🛠️ For questions or issues, contact ",
                       a("John", href = "https://github.com/jcaperella29", target = "_blank"),
                       " or view the full documentation at ",
                       a("GWAS_SHINY_APP", href = "https://github.com/jcaperella29/GWAS_SHINY_APP", target = "_blank"), "."
                     )
                   )
                 )
        ),
        
        
        tabPanel("Phenotype Preview", DTOutput("phenoPreview")),
        tabPanel("Covariates Preview", DTOutput("covarPreview")),
        tabPanel("QC Summary", DTOutput("qcTable")),
        tabPanel("Filtered SNPs", DTOutput("filteredSnps")),
        tabPanel("PCA Plot", plotlyOutput("pcaPlot")),
        
        tabPanel("UMAP Plot", plotlyOutput("umapPlot")),
        
        tabPanel("Power Analysis", DTOutput("powerTable")),
        
        tabPanel("Power Curve", plotlyOutput("powerPlot")),
        
        tabPanel("GWAS Results", DTOutput("gwasResults")),
        tabPanel("Covariate Effects", 
                 DTOutput("covarPvals")),
        tabPanel("Final Phenotype Hits", 
                 DTOutput("finalHits")),
        tabPanel("Q-Q Plot", plotlyOutput("qqPlot")),
        tabPanel("Manhattan Plot", plotlyOutput("manhattanPlot")),
        tabPanel("SNP-to_Gene Mapping",
                 DTOutput("annotatedHits")),
        tabPanel("Gene Map Plot",
                 plotlyOutput("geneMapPlot", height = "500px")),
        tabPanel("Pathway Enrichment",
                 DTOutput("enrichmentTable")),
        
        tabPanel("Enrichment Barplot", plotlyOutput("enrichBarplot")),
        tabPanel("Random Forest Results",
                 tabsetPanel(
                   tabPanel("Predictions", DTOutput("rfPredictions")),
                   tabPanel("Metrics", DTOutput("rfMetrics")),
                   tabPanel("Variable Importance", DTOutput("rfImportance"))
                 )
        ),
        tabPanel("ROC Curve", plotlyOutput("rfROC"))
        
        
      )
      
      
    )
  )
)


server <- function(input, output, session) {
  
  rv <- reactiveValues(
    geno_df = NULL,
    pheno = NULL,
    covar = NULL,
    qc_snps = NULL,
    qc_table = NULL,
    gwas_res = NULL,
    power_df = NULL,
    
    ml_out = NULL
  )
  
  
  output$qqPlot <- renderPlotly(NULL)
  output$manhattanPlot <- renderPlotly(NULL)
  
  observeEvent(input$phenoFile, {
    rv$pheno <- fread(input$phenoFile$datapath)
    output$phenoPreview <- renderDT(rv$pheno)
  })
  
  observeEvent(input$covarFile, {
    rv$covar <- fread(input$covarFile$datapath)
    output$covarPreview <- renderDT(rv$covar)
  })
  
  
  observeEvent(input$vcfFile, {
    vcf <- readVcf(input$vcfFile$datapath, genome = "hg19")
    gt <- geno(vcf)$GT
    
    dosage <- apply(gt, 2, function(x) {
      sapply(x, function(g) {
        if (g %in% c("0/0", "0|0")) return(0)
        if (g %in% c("0/1", "1/0", "0|1", "1|0")) return(1)
        if (g %in% c("1/1", "1|1")) return(2)
        return(NA_real_)
      })
    })
    
    snp_ids <- make.names(names(rowRanges(vcf)), unique = TRUE)
    rownames(dosage) <- snp_ids
    colnames(dosage) <- colnames(gt)
    
    # Store metadata vectors indexed by SNP ID
    snp_chr_vec <- as.character(seqnames(rowRanges(vcf)))
    snp_pos_vec <- start(rowRanges(vcf))
    
    names(snp_chr_vec) <- snp_ids
    names(snp_pos_vec) <- snp_ids
    
    rv$geno_df <- as.data.frame(dosage)
    rv$snp_chr <- snp_chr_vec  # named vector: SNP -> CHR
    rv$snp_pos <- snp_pos_vec  # named vector: SNP -> POS
    
    
    # After this:
    rv$snp_chr <- snp_chr_vec
    rv$snp_pos <- snp_pos_vec
    
    # Add this:
    rv$snp_info <- data.table(
      SNP = names(snp_chr_vec),
      CHR = snp_chr_vec,
      POS = snp_pos_vec
    )
    
  })
  
  
  
  
  
  observeEvent(input$runQC, {
    req(rv$geno_df)
    showNotification("Running QC...", type = "message", duration = NULL, id = "qc_notice")
    
    geno <- rv$geno_df
    snp_ids <- rownames(geno)
    maf <- rowMeans(geno, na.rm = TRUE) / 2
    call_rate <- rowMeans(!is.na(geno))
    allele_freq <- maf * 2
    
    hwe_pvals <- sapply(1:nrow(geno), function(i) {
      g <- geno[i, ]
      g <- g[!is.na(g)]
      obs <- table(factor(g, levels = c(0, 1, 2)))
      n <- sum(obs)
      if (n == 0) return(NA_real_)
      p <- mean(g) / 2
      q <- 1 - p
      exp0 <- n * q^2
      exp1 <- n * 2 * p * q
      exp2 <- n * p^2
      expected <- c(exp0, exp1, exp2)
      if (any(expected < 1)) return(NA_real_)
      chisq <- sum((obs - expected)^2 / expected)
      pchisq(chisq, df = 1, lower.tail = FALSE)
    })
    
    pass <- maf > input$maf_thresh &
      allele_freq > input$af_thresh &
      call_rate > 0.95 &
      hwe_pvals > input$hwe_thresh
    pass[is.na(pass)] <- FALSE
    
    qc_table <- data.table(
      SNP = snp_ids,
      MAF = round(maf, 4),
      AF = round(allele_freq, 4),
      CallRate = round(call_rate, 4),
      HWE_P = signif(hwe_pvals, 4),
      Pass = pass
    )
    
    rv$qc_table <- qc_table
    rv$qc_snps <- geno[qc_table$Pass, , drop = FALSE]
    
    
    output$qcTable <- renderDT(qc_table)
    output$filteredSnps <- renderDT({
      if (!is.null(rv$qc_snps)) {
        data.table(SNP = rownames(rv$qc_snps))
      }
    })
    
    removeNotification("qc_notice")
    showNotification(paste0("QC complete. Passed SNPs: ", sum(pass, na.rm = TRUE)), type = "message", duration = 5)
  })
  
  
  observeEvent(input$runGWAS, {
    req(rv$qc_snps, rv$pheno)
    
    showNotification("🚀 Running GWAS...", type = "message", duration = NULL, id = "gwas_notice")
    
    geno_df <- as.data.frame(t(rv$qc_snps))  # samples x SNPs
    pheno <- rv$pheno
    covar <- rv$covar
    
    # Align sample IDs
    common_ids <- rownames(geno_df)
    if (!is.null(covar)) {
      common_ids <- Reduce(intersect, list(rownames(geno_df), pheno$sample_id, covar$sample_id))
      covar <- covar[sample_id %in% common_ids]
      covar <- covar[match(common_ids, covar$sample_id), ]
    } else {
      common_ids <- intersect(rownames(geno_df), pheno$sample_id)
    }
    
    if (length(common_ids) == 0) {
      showNotification("❌ No matching sample IDs found!", type = "error")
      removeNotification("gwas_notice")
      return()
    }
    
    geno_df <- geno_df[common_ids, , drop = FALSE]
    pheno <- pheno[match(common_ids, pheno$sample_id), ]
    
    # GWAS loop
    gwas_res <- data.table()
    na_skipped <- 0
    mono_skipped <- 0
    model_failed <- 0
    
    for (snp in colnames(geno_df)) {
      g <- geno_df[[snp]]
      if (all(is.na(g))) { na_skipped <- na_skipped + 1; next }
      if (length(unique(g[!is.na(g)])) <= 1) { mono_skipped <- mono_skipped + 1; next }
      
      model <- tryCatch(glm(pheno$trait ~ g, family = binomial()), error = function(e) NULL)
      if (!is.null(model)) {
        pval <- tryCatch(summary(model)$coefficients[2, 4], error = function(e) NA)
        gwas_res <- rbind(gwas_res, data.table(
          SNP = snp,
          CHR = rv$snp_chr[[snp]],
          POS = rv$snp_pos[[snp]],
          P = pval
        ))
      } else {
        model_failed <- model_failed + 1
      }
    }
    
    if (nrow(gwas_res) == 0) {
      showNotification("❌ GWAS complete but no SNPs modeled!", type = "error")
      showNotification(paste("⛔ NA:", na_skipped, "| Mono:", mono_skipped, "| Fail:", model_failed), type = "error")
      removeNotification("gwas_notice")
      return()
    }
    
    # 💥 Bonferroni correction
    m_tests <- nrow(gwas_res)
    gwas_res[, P_bonf := pmin(P * m_tests, 1)]
    
    rv$gwas_res <- gwas_res
    
    output$gwasResults <- renderDT({
      req(rv$gwas_res)
      
      # Choose filter column
      filter_col <- if (input$use_bonferroni) "P_bonf" else "P"
      pvals <- rv$gwas_res[[filter_col]]
      
      # Apply threshold
      filtered <- rv$gwas_res[!is.na(pvals) & pvals < input$pval_thresh]
      
      if (nrow(filtered) == 0) {
        showNotification(paste("⚠️ No SNPs passed", if (input$use_bonferroni) "Bonferroni" else "raw", "threshold."), type = "warning")
        return(rv$gwas_res[order(P)])
      }
      
      filtered[order(get(filter_col))]
    })
    
    
    # Covariate association testing
    covar_assoc <- data.table()
    if (!is.null(covar)) {
      for (snp in colnames(geno_df)) {
        g <- geno_df[[snp]]
        if (length(unique(g[!is.na(g)])) <= 1) next
        if (any(is.na(g))) next
        
        for (cov in setdiff(names(covar), "sample_id")) {
          cv <- covar[[cov]]
          model <- tryCatch(glm(cv ~ g, family = gaussian()), error = function(e) NULL)
          if (!is.null(model)) {
            pval <- tryCatch(summary(model)$coefficients[2, 4], error = function(e) NA)
            if (!is.na(pval)) {
              covar_assoc <- rbind(covar_assoc, data.table(
                SNP = snp,
                CHR = rv$snp_chr[[snp]],
                POS = rv$snp_pos[[snp]],
                Covariate = cov,
                P_value = signif(pval, 4)
              ))
            }
          }
        }
      }
    }
    
    rv$covar_assoc <- covar_assoc
    output$covarPvals <- renderDT(covar_assoc)
    
    # Final hits (phenotype-only)
    sig_snps <- gwas_res[P_bonf < input$pval_thresh]$SNP
    covar_snps <- unique(covar_assoc[P_value < input$pval_thresh]$SNP)
    final_snps <- setdiff(sig_snps, covar_snps)
    
    final_hits <- gwas_res[SNP %in% final_snps]
    setnames(final_hits, "P", "P_pheno")
    final_hits <- final_hits[, .(SNP, CHR, POS, P_pheno, P_bonf)]
    
    rv$final_hits <- final_hits
    output$finalHits <- renderDT({
      req(rv$final_hits)
      validate(need(nrow(rv$final_hits) > 0, "⚠️ No phenotype-only significant SNPs"))
      rv$final_hits
    })
    
    removeNotification("gwas_notice")
    showNotification(paste0("✅ GWAS done. SNPs: ", nrow(gwas_res), " | NA:", na_skipped,
                            " | Mono:", mono_skipped, " | Fail:", model_failed), type = "message", duration = 5)
  })
  
  
  
  
  
  observeEvent(input$runQQ, {
    req(rv$gwas_res)
    showNotification("Generating Q-Q plot...", type = "message", duration = NULL, id = "qq_notice")
    
    tryCatch({
      pvals <- rv$gwas_res$P
      pvals <- pvals[!is.na(pvals) & pvals > 0]
      
      if (length(pvals) < 3) {
        output$qqPlot <- renderPlotly({
          plot_ly() %>% layout(title = "Q-Q Plot",
                               annotations = list(text = "Not enough SNPs to plot", showarrow = FALSE))
        })
        removeNotification("qq_notice")
        return(showNotification("Q-Q plot skipped: too few SNPs", type = "warning", duration = 4))
      }
      
      obs <- -log10(sort(pvals))
      exp <- -log10(ppoints(length(obs)))
      
      output$qqPlot <- renderPlotly({
        plot_ly(x = exp, y = obs, type = "scatter", mode = "markers") %>%
          layout(title = "Q-Q Plot",
                 xaxis = list(title = "Expected -log10(p)"),
                 yaxis = list(title = "Observed -log10(p)"),
                 shapes = list(
                   list(type = "line", x0 = min(exp), x1 = max(exp),
                        y0 = min(exp), y1 = max(exp),
                        line = list(dash = "dash", color = "gray"))
                 ))
      })
      removeNotification("qq_notice")
      showNotification("Q-Q plot complete.", type = "message", duration = 4)
    }, error = function(e) {
      removeNotification("qq_notice")
      showNotification(paste("Q-Q plot failed:", e$message), type = "error", duration = 6)
    })
  })
  
  observeEvent(input$runManhattan, {
    req(rv$gwas_res)
    
    # Ensure snp_info is present (or reconstruct it from SNP metadata)
    if (is.null(rv$snp_info) && !is.null(rv$snp_chr) && !is.null(rv$snp_pos)) {
      rv$snp_info <- data.table(
        SNP = names(rv$snp_chr),
        CHR = rv$snp_chr,
        POS = rv$snp_pos
      )
    }
    
    req(rv$snp_info)
    
    showNotification("🔬 Generating Manhattan plot...", type = "message", duration = NULL, id = "manhattan_notice")
    
    tryCatch({
      # Merge GWAS results with SNP info (CHR + POS from VCF)
      gwas <- merge(rv$gwas_res, rv$snp_info, by = "SNP", all.x = TRUE)
      
      # Resolve possible merge duplicates
      if ("CHR.x" %in% names(gwas)) gwas[, CHR := CHR.x]
      if ("POS.x" %in% names(gwas)) gwas[, POS := POS.x]
      gwas[, c("CHR.x", "CHR.y", "POS.x", "POS.y") := NULL]
      
      # Filter and validate
      gwas <- gwas[!is.na(P) & P > 0 & !is.na(CHR) & !is.na(POS)]
      if (nrow(gwas) < 3) {
        output$manhattanPlot <- renderPlotly({
          plot_ly() %>% layout(
            title = "Manhattan Plot",
            annotations = list(text = "Not enough SNPs to plot", showarrow = FALSE)
          )
        })
        removeNotification("manhattan_notice")
        return(showNotification("Manhattan plot skipped: too few SNPs", type = "warning", duration = 4))
      }
      
      # Format types
      gwas[, CHR := as.character(CHR)]
      gwas[, POS := as.numeric(POS)]
      gwas[, logP := -log10(P)]
      
      # Region filter
      if (input$chr_filter != "All") {
        gwas <- gwas[CHR == input$chr_filter]
      }
      gwas <- gwas[POS >= input$pos_min & POS <= input$pos_max]
      
      if (nrow(gwas) < 3) {
        output$manhattanPlot <- renderPlotly({
          plot_ly() %>% layout(
            title = "Manhattan Plot",
            annotations = list(text = "No SNPs in selected region", showarrow = FALSE)
          )
        })
        removeNotification("manhattan_notice")
        return(showNotification("No SNPs in selected region", type = "warning", duration = 4))
      }
      
      # Plot
      sig_thresh <- -log10(input$pval_thresh)
      
      output$manhattanPlot <- renderPlotly({
        plot_ly(
          gwas,
          x = ~POS,
          y = ~logP,
          type = "scatter",
          mode = "markers",
          marker = list(color = "#1f77b4", size = 6),
          text = ~paste0(
            "SNP: ", SNP, "<br>",
            "Chr: ", CHR, "<br>",
            "Pos: ", POS, "<br>",
            "P: ", formatC(P, format = "e", digits = 2)
          ),
          hoverinfo = "text"
        ) %>%
          layout(
            title = paste("Manhattan Plot - Chr", ifelse(input$chr_filter == "All", "All", input$chr_filter)),
            xaxis = list(title = "Genomic Position"),
            yaxis = list(title = "-log10(p)"),
            shapes = list(
              list(
                type = "line",
                x0 = min(gwas$POS),
                x1 = max(gwas$POS),
                y0 = sig_thresh,
                y1 = sig_thresh,
                line = list(color = "red", dash = "dash")
              )
            )
          )
      })
      
      removeNotification("manhattan_notice")
      showNotification("✅ Manhattan plot complete.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("manhattan_notice")
      showNotification(paste("❌ Manhattan plot failed:", e$message), type = "error", duration = 6)
    })
  })
  
  
  
  observeEvent(input$runML, {
    req(rv$gwas_res, rv$pheno, rv$qc_snps)
    sig_snps <- rv$gwas_res[P < 5e-8]$SNP
    if (length(sig_snps) < 1) {
      output$mlResults <- renderDT(data.table(Message = "No significant SNPs"))
      return()
    }
    
    geno_df <- rv$qc_snps[sig_snps, , drop = FALSE]
    geno_df <- as.data.frame(t(geno_df))
    pheno <- rv$pheno
    common_ids <- intersect(rownames(geno_df), pheno$sample_id)
    geno_df <- geno_df[common_ids, , drop = FALSE]
    pheno <- pheno[match(common_ids, pheno$sample_id), ]
    ml_data <- data.table(geno_df, trait = as.factor(pheno$trait))
    
    task <- TaskClassif$new("gwas_class", backend = ml_data, target = "trait")
    learner <- lrn("classif.ranger", importance = "impurity")
    learner$train(task)
    imp <- learner$importance()
    rv$ml_out <- data.table(SNP = names(imp), Importance = imp)
    output$mlResults <- renderDT(rv$ml_out)
  })
  
  # GWAS full results
  output$download_gwas <- downloadHandler(
    filename = function() {
      paste0("GWAS_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fwrite(rv$gwas_res, file)
    }
  )
  
  # Covariate associations
  output$download_covar <- downloadHandler(
    filename = function() {
      paste0("Covariate_Pvals_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fwrite(rv$covar_assoc, file)
    }
  )
  
  # Final hits (phenotype-only SNPs)
  output$download_final <- downloadHandler(
    filename = function() {
      paste0("Final_Hits_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fwrite(rv$final_hits, file)
    }
  )
  
  observeEvent(input$mapGenes, {
    req(rv$gwas_res)
    
    showNotification("🔍 Mapping SNPs to nearest genes...", type = "message", duration = NULL, id = "gene_notice")
    
    tryCatch({
      library(GenomicFeatures)
      library(org.Hs.eg.db)
      library(TxDb.Hsapiens.UCSC.hg19.knownGene)
      library(AnnotationDbi)
      
      # Determine which column to filter on
      filter_col <- if (input$use_bonferroni) "P_bonf" else "P"
      pvals <- rv$gwas_res[[filter_col]]
      
      # Apply threshold to define SNPs to map
      hits <- rv$gwas_res[!is.na(pvals) & pvals < input$pval_thresh]
      
      if (nrow(hits) < 1) {
        removeNotification("gene_notice")
        showNotification("⚠️ No SNPs passed threshold for gene mapping.", type = "warning", duration = 5)
        return()
      }
      
      # SNP GRanges
      snp_gr <- GRanges(
        seqnames = hits$CHR,
        ranges = IRanges(start = hits$POS, width = 1),
        SNP = hits$SNP
      )
      
      genes_gr <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene, single.strand.genes.only = FALSE)
      genes_gr <- unlist(genes_gr)
      
      # Match styles
      seqlevelsStyle(snp_gr) <- seqlevelsStyle(genes_gr)
      
      # Nearest genes
      nearest_idxs <- nearest(snp_gr, genes_gr)
      matched_genes <- genes_gr[nearest_idxs]
      
      # Gene IDs + symbols
      entrez_ids <- names(matched_genes)
      gene_symbols <- mapIds(org.Hs.eg.db,
                             keys = entrez_ids,
                             column = "SYMBOL",
                             keytype = "ENTREZID",
                             multiVals = "first")
      
      # Distance calc
      dist_to_gene <- abs(start(matched_genes) - start(snp_gr))
      
      # Annotated result
      annotated_hits <- data.table(
        hits,
        Gene_ID = entrez_ids,
        Gene_Symbol = gene_symbols[entrez_ids],
        Gene_Distance = dist_to_gene
      )
      
      rv$annotated_hits <- annotated_hits
      
      output$annotatedHits <- renderDT({
        req(rv$annotated_hits)
        annotated_hits
      })
      
      output$download_annotated <- downloadHandler(
        filename = function() paste0("Annotated_Hits_", Sys.Date(), ".csv"),
        content = function(file) {
          fwrite(rv$annotated_hits, file)
        }
      )
      
      removeNotification("gene_notice")
      showNotification(paste0("✅ Mapped ", nrow(annotated_hits), " SNPs to nearest genes."), type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("gene_notice")
      showNotification(paste("❌ Gene mapping failed:", e$message), type = "error", duration = 6)
    })
  })
  
  
  observeEvent(input$plotGenes, {
    req(rv$annotated_hits)
    
    showNotification("🧬 Generating gene map...", type = "message", duration = NULL, id = "gene_map_notice")
    
    tryCatch({
      hits <- rv$annotated_hits
      
      output$geneMapPlot <- renderPlotly({
        plot_ly(hits, x = ~POS, y = ~Gene_Distance, type = 'scatter', mode = 'markers+text',
                text = ~paste0(
                  "<b>SNP:</b> ", SNP, "<br>",
                  "<b>Gene:</b> ", Gene_Symbol, "<br>",
                  "<b>Chr:</b> ", CHR, "<br>",
                  "<b>Pos:</b> ", POS, "<br>",
                  "<b>Distance:</b> ", Gene_Distance, " bp"
                ),
                hoverinfo = 'text',
                textposition = "top center",
                marker = list(size = 10, color = 'blue')) %>%
          layout(
            title = "Gene Mapping: SNP → Nearest Gene",
            xaxis = list(title = "SNP Position (bp)"),
            yaxis = list(title = "Distance to Nearest Gene (bp)"),
            margin = list(t = 70)
          )
      })
      
      removeNotification("gene_map_notice")
      showNotification("✅ Gene map ready.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("gene_map_notice")
      showNotification(paste("❌ Gene plot failed:", e$message), type = "error", duration = 6)
    })
  })
  
  observeEvent(input$runEnrichment, {
    req(rv$annotated_hits)
    
    showNotification("📡 Running gene set enrichment...", type = "message", duration = NULL, id = "enrich_notice")
    
    tryCatch({
      genes <- unique(na.omit(rv$annotated_hits$Gene_Symbol))
      if (length(genes) < 2) stop("Need at least 2 gene symbols to enrich.")
      
      enrichr_dbs <- c("KEGG_2021_Human", "GO_Biological_Process_2021", "Reactome_2016")
      enrichr_results <- enrichr(genes, enrichr_dbs)
      
      # Combine all DBs into one big table
      enrichment_combined <- rbindlist(
        lapply(names(enrichr_results), function(db) {
          dt <- as.data.table(enrichr_results[[db]])
          dt[, Database := db]
          return(dt)
        }),
        fill = TRUE
      )
      
      rv$enrichment_results <- enrichment_combined
      
      output$enrichmentTable <- renderDT({
        datatable(enrichment_combined, options = list(pageLength = 10, scrollX = TRUE))
      })
      
      output$download_enrichment <- downloadHandler(
        filename = function() paste0("Pathway_Enrichment_", Sys.Date(), ".csv"),
        content = function(file) {
          fwrite(rv$enrichment_results, file)
        }
      )
      
      removeNotification("enrich_notice")
      showNotification("✅ Enrichment complete.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("enrich_notice")
      showNotification(paste("❌ Enrichment failed:", e$message), type = "error", duration = 6)
    })
  })
  
  observeEvent(input$plotEnrichment, {
    req(rv$enrichment_results)
    
    showNotification("📊 Plotting enrichment results...", type = "message", duration = NULL, id = "plot_notice")
    
    tryCatch({
      selected_db <- input$enrich_db
      enrich_df <- rv$enrichment_results[Database == selected_db]
      
      if (nrow(enrich_df) == 0) {
        removeNotification("plot_notice")
        showNotification(paste("⚠️ No enrichment results for", selected_db), type = "warning")
        return()
      }
      
      enrich_df <- enrich_df[order(Adjusted.P.value)]
      top_terms <- head(enrich_df, 10)
      
      output$enrichBarplot <- renderPlotly({
        plot_ly(
          data = top_terms,
          x = ~reorder(as.character(Term), -as.numeric(Adjusted.P.value)),
          y = ~-log10(as.numeric(Adjusted.P.value)),
          type = 'bar',
          text = ~paste("Genes:", Genes),
          hoverinfo = 'text+y',
          marker = list(color = 'darkgreen')
        ) %>%
          layout(
            title = paste("Top Enriched Pathways:", selected_db),
            xaxis = list(title = "Pathway", tickangle = -45),
            yaxis = list(title = "-log10(Adj P-value)"),
            margin = list(b = 100)
          )
      })
      
      removeNotification("plot_notice")
      showNotification("✅ Barplot generated.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("plot_notice")
      showNotification(paste("❌ Plot error:", e$message), type = "error", duration = 6)
    })
  })
  
  observeEvent(input$runPCA, {
    req(rv$qc_snps)
    
    showNotification("📈 Running PCA on QC-passed SNPs...", type = "message", duration = NULL, id = "pca_notice")
    
    tryCatch({
      pca_data <- t(rv$qc_snps)
      pca_res <- prcomp(pca_data, center = TRUE, scale. = TRUE)
      
      pcs <- as.data.table(pca_res$x[, 1:2])
      pcs[, Sample := rownames(pca_res$x)]
      
      if (!is.null(rv$pheno)) {
        pcs <- merge(pcs, rv$pheno[, .(sample_id, trait)], by.x = "Sample", by.y = "sample_id", all.x = TRUE)
      }
      
      output$pcaPlot <- renderPlotly({
        plot_ly(
          data = pcs,
          x = ~PC1,
          y = ~PC2,
          color = ~trait,
          text = ~Sample,
          type = "scatter",
          mode = "markers",
          marker = list(size = 10)
        ) %>%
          layout(
            title = "PCA of Genotype Matrix",
            xaxis = list(title = "PC1"),
            yaxis = list(title = "PC2"),
            hovermode = "closest"
          )
      })
      
      removeNotification("pca_notice")
      showNotification("✅ PCA complete.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("pca_notice")
      showNotification(paste("❌ PCA failed:", e$message), type = "error", duration = 6)
    })
  })
  observeEvent(input$runUMAP, {
    req(rv$qc_snps)
    
    showNotification("📉 Running UMAP on QC-passed SNPs...", type = "message", duration = NULL, id = "umap_notice")
    
    tryCatch({
      umap_data <- t(rv$qc_snps)
      
      # Auto-install UMAP if needed
      if (!requireNamespace("umap", quietly = TRUE)) {
        install.packages("umap")
      }
      library(umap)
      
      umap_res <- umap(umap_data)
      
      umap_dt <- as.data.table(umap_res$layout)
      setnames(umap_dt, c("UMAP1", "UMAP2"))
      umap_dt[, Sample := rownames(umap_data)]
      
      if (!is.null(rv$pheno)) {
        umap_dt <- merge(umap_dt, rv$pheno[, .(sample_id, trait)], by.x = "Sample", by.y = "sample_id", all.x = TRUE)
      }
      
      output$umapPlot <- renderPlotly({
        plot_ly(
          data = umap_dt,
          x = ~UMAP1,
          y = ~UMAP2,
          color = ~trait,
          text = ~Sample,
          type = "scatter",
          mode = "markers",
          marker = list(size = 10)
        ) %>%
          layout(
            title = "UMAP of Genotype Matrix",
            xaxis = list(title = "UMAP1"),
            yaxis = list(title = "UMAP2"),
            hovermode = "closest"
          )
      })
      
      removeNotification("umap_notice")
      showNotification("✅ UMAP complete.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("umap_notice")
      showNotification(paste("❌ UMAP failed:", e$message), type = "error", duration = 6)
    })
  })
  
  observeEvent(input$runPower, {
    req(rv$pheno, rv$qc_snps)
    
    showNotification("🧮 Running power analysis...", type = "message", duration = NULL, id = "power_notice")
    
    tryCatch({
      library(pwr)
      
      geno <- as.data.frame(t(rv$qc_snps))  # samples x SNPs
      common_ids <- intersect(rownames(geno), rv$pheno$sample_id)
      
      geno <- geno[common_ids, , drop = FALSE]
      pheno <- rv$pheno[match(common_ids, rv$pheno$sample_id), ]
      geno$trait <- pheno$trait
      
      if (!is.numeric(geno$trait) || !all(geno$trait %in% c(0, 1))) {
        removeNotification("power_notice")
        showNotification("⚠️ Trait must be binary (0/1) for power analysis.", type = "warning")
        return()
      }
      
      x0 <- geno[geno$trait == 0, 1]
      x1 <- geno[geno$trait == 1, 1]
      
      if (length(x0) < 2 || length(x1) < 2) {
        removeNotification("power_notice")
        showNotification("⚠️ Not enough samples in one or both groups.", type = "warning")
        return()
      }
      
      # Cohen's d
      d <- abs(mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) /
        sqrt((var(x0, na.rm = TRUE) + var(x1, na.rm = TRUE)) / 2)
      
      # Observed power
      pw <- pwr.t.test(
        n = min(length(x0), length(x1)),
        d = d,
        sig.level = 0.05,
        type = "two.sample",
        alternative = "two.sided"
      )
      
      # Required N per group for 80% power
      required_n <- pwr.t.test(
        power = 0.8,
        d = d,
        sig.level = 0.05,
        type = "two.sample",
        alternative = "two.sided"
      )$n
      
      power_df <- data.table(
        SNP = colnames(geno)[1],
        Group0_N = length(x0),
        Group1_N = length(x1),
        Cohen_d = round(d, 3),
        Observed_Power = round(pw$power, 3),
        Required_N_80Power = ceiling(required_n)
      )
      rv$power_df <- power_df
      
      output$powerTable <- renderDT({
        datatable(rv$power_df, options = list(dom = 't', scrollX = TRUE))
      })
      
      removeNotification("power_notice")
      showNotification("✅ Power analysis complete.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("power_notice")
      showNotification(paste("❌ Power analysis failed:", e$message), type = "error", duration = 6)
    })
  })
  
      
   
  
  output$download_power <- downloadHandler(
    filename = function() {
      paste0("Power_Analysis_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(rv$power_df)
      fwrite(rv$power_df, file)
    }
  )
  
  
  observeEvent(input$runPowerCurve, {
    req(rv$power_df)
    
    showNotification("📈 Generating power curve...", type = "message", duration = NULL, id = "powercurve_notice")
    
    tryCatch({
      d <- rv$power_df$Cohen_d[1]
      if (is.null(d) || is.na(d) || d <= 0) {
        removeNotification("powercurve_notice")
        showNotification("⚠️ Invalid effect size for power curve.", type = "warning")
        return()
      }
      
      n_range <- seq(5, 200, by = 1)
      power_vals <- sapply(n_range, function(n) {
        tryCatch({
          pwr.t.test(n = n, d = d, sig.level = 0.05,
                     type = "two.sample", alternative = "two.sided")$power
        }, error = function(e) NA)
      })
      
      power_df <- data.table(SampleSizePerGroup = n_range, Power = power_vals)
      
      output$powerPlot <- renderPlotly({
        plot_ly(power_df, x = ~SampleSizePerGroup, y = ~Power,
                type = 'scatter', mode = 'lines', line = list(color = "blue")) %>%
          layout(
            title = "Power Curve (based on Cohen's d)",
            xaxis = list(title = "Sample Size Per Group"),
            yaxis = list(title = "Statistical Power", range = c(0, 1)),
            shapes = list(
              list(type = "line", x0 = min(n_range), x1 = max(n_range),
                   y0 = 0.8, y1 = 0.8,
                   line = list(dash = "dash", color = "red"))
            ),
            annotations = list(
              list(x = max(n_range) * 0.85, y = 0.82,
                   text = "80% Power Threshold", showarrow = FALSE,
                   font = list(color = "red", size = 12))
            )
          )
      })
      
      removeNotification("powercurve_notice")
      showNotification("✅ Power curve ready.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("powercurve_notice")
      showNotification(paste("❌ Power curve failed:", e$message), type = "error", duration = 6)
    })
  })
  observeEvent(input$runRF, {
    req(rv$final_hits, rv$qc_snps, rv$pheno)
    
    showNotification("🌲 Training Random Forest model...", type = "message", duration = NULL, id = "rf_notice")
    
    tryCatch({
      final_snps <- rv$final_hits$SNP
      if (length(final_snps) < 2) {
        removeNotification("rf_notice")
        showNotification("⚠️ Not enough SNPs in final hits for ML.", type = "warning")
        return()
      }
      
      geno <- as.data.frame(t(rv$qc_snps[final_snps, , drop = FALSE]))
      common_ids <- intersect(rownames(geno), rv$pheno$sample_id)
      geno <- geno[common_ids, , drop = FALSE]
      pheno <- rv$pheno[match(common_ids, rv$pheno$sample_id), ]
      
      # Ensure binary trait
      if (!all(pheno$trait %in% c(0, 1))) {
        removeNotification("rf_notice")
        showNotification("⚠️ Trait must be binary (0/1).", type = "warning")
        return()
      }
      
      df <- data.table(geno, trait = as.factor(pheno$trait))
      
      # Split train/test 70/30
      set.seed(123)
      train_idx <- sample(seq_len(nrow(df)), size = 0.7 * nrow(df))
      train <- df[train_idx]
      test <- df[-train_idx]
      
      task <- TaskClassif$new(id = "rf_task", backend = train, target = "trait")
      learner <- lrn("classif.ranger", predict_type = "prob", importance = "impurity")
      learner$train(task)
      
      # Predict
      pred <- learner$predict_newdata(test)
      
      # Predictions table
      preds_df <- data.table(
        Sample = rownames(test),
        True = test$trait,
        Pred = pred$response,
        Prob_1 = round(pred$prob[, "1"], 4)
      )
      
      output$rfPredictions <- renderDT({
        datatable(preds_df, options = list(pageLength = 10))
      })
      rv$rf_preds <- preds_df
      
      # Metrics
      cm <- table(True = preds_df$True, Pred = preds_df$Pred)
      TP <- cm["1", "1"]
      TN <- cm["0", "0"]
      FP <- cm["0", "1"]
      FN <- cm["1", "0"]
      
      acc <- (TP + TN) / sum(cm)
      sens <- TP / (TP + FN)
      spec <- TN / (TN + FP)
      
      auc <- tryCatch({
        pred_obj <- ROCR::prediction(preds_df$Prob_1, preds_df$True)
        ROCR::performance(pred_obj, "auc")@y.values[[1]]
      }, error = function(e) NA)
      
      metrics_df <- data.table(
        Metric = c("Accuracy", "Sensitivity", "Specificity", "AUC"),
        Value = round(c(acc, sens, spec, auc), 4)
      )
      
      output$rfMetrics <- renderDT({
        datatable(metrics_df)
      })
      
      # Variable Importance
      vi <- learner$importance()
      vi_df <- data.table(SNP = names(vi), Importance = round(vi, 4))
      vi_df <- vi_df[order(-Importance)]
      
      output$rfImportance <- renderDT({
        datatable(vi_df)
      })
      
      
      
      showNotification("✅ Random Forest model trained.", type = "message", duration = 4)
      removeNotification("rf_notice")
      
    }, error = function(e) {
      removeNotification("rf_notice")
      showNotification(paste("❌ RF training failed:", e$message), type = "error", duration = 6)
    })
  })
  observeEvent(input$runROC, {
    req(rv$rf_preds, rv$final_hits, rv$qc_snps, rv$pheno)
    if (length(rv$final_hits$SNP) < 2) {
      showNotification("⚠️ Not enough SNPs in final hits to build ROC.", type = "warning")
      return()
    }
    
    showNotification("📉 Generating ROC plot...", type = "message", duration = NULL, id = "roc_notice")
    
    tryCatch({
      # Reconstruct prediction from stored df
      truth <- as.numeric(as.character(rv$rf_preds$True))
      prob_1 <- rv$rf_preds$Prob_1
      
      pred_obj <- ROCR::prediction(prob_1, truth)
      perf <- ROCR::performance(pred_obj, "tpr", "fpr")
      
      fpr <- perf@x.values[[1]]
      tpr <- perf@y.values[[1]]
      roc_df <- data.table(FPR = fpr, TPR = tpr)
      
      output$rfROC <- renderPlotly({
        plot_ly(roc_df, x = ~FPR, y = ~TPR, type = 'scatter', mode = 'lines',
                line = list(color = "darkgreen")) %>%
          layout(
            title = "ROC Curve",
            xaxis = list(title = "False Positive Rate"),
            yaxis = list(title = "True Positive Rate"),
            shapes = list(list(
              type = "line", x0 = 0, x1 = 1, y0 = 0, y1 = 1,
              line = list(dash = "dash", color = "gray")
            ))
          )
      })
      
      removeNotification("roc_notice")
      showNotification("✅ ROC curve ready.", type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("roc_notice")
      showNotification(paste("❌ ROC plot failed:", e$message), type = "error", duration = 6)
    })
  })
  output$download_rf_preds <- downloadHandler(
    filename = function() paste0("RF_Predictions_", Sys.Date(), ".csv"),
    content = function(file) {
      req(rv$rf_preds)
      fwrite(rv$rf_preds, file)
    }
  )
  
  output$download_rf_metrics <- downloadHandler(
    filename = function() paste0("RF_Metrics_", Sys.Date(), ".csv"),
    content = function(file) {
      req(output$rfMetrics)
      fwrite(as.data.table(output$rfMetrics()), file)
    }
  )
  
  output$download_rf_importance <- downloadHandler(
    filename = function() paste0("RF_Importance_", Sys.Date(), ".csv"),
    content = function(file) {
      req(output$rfImportance)
      fwrite(as.data.table(output$rfImportance()), file)
    }
  )
  
  
  
}



shinyApp(ui, server)
