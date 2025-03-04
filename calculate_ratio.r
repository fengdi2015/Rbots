# Load necessary library
library(dplyr)

# Define gene symbols
genes <- c("TP53", "BRCA1", "EGFR", "KRAS", "PTEN")

# Generate mock data
set.seed(123)  # For reproducibility
data <- data.frame(
  Sample_ID = rep(paste0("Sample_", 1:20), each = sample(1:3, 20, replace = TRUE)),
  Gene_Symbol = sample(genes, 50, replace = TRUE),
  Mutation_Type = sample(c("S", "NS"), 50, replace = TRUE)
)

# Count occurrences of S and NS per Sample ID and Gene Symbol
mutation_counts <- data %>%
  group_by(Sample_ID, Gene_Symbol, Mutation_Type) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Mutation_Type, values_from = Count, values_fill = list(Count = 0))

# Ensure both S and NS columns exist
mutation_counts <- mutation_counts %>%
  mutate(S = ifelse(is.na(S), 0, S),
         NS = ifelse(is.na(NS), 0, NS),
         S_NS_Ratio = S / (NS + 1e-6))  # Avoid division by zero

# Print result
print(mutation_counts)
