# 🧬 DNA Folding Simulator

A modular simulation that generates plausible DNA strand structures using trend-following stochastic modeling. This tool creates semi-random nucleotide positions that follow a general spatial trajectory, adding local noise to reflect real-world irregularity in folding behavior.

---

## 📖 Overview
This simulator approximates DNA folding using simplified physical principles. It does not perform full molecular dynamics, but incorporates the following key interactions:

- **Backbone Bonds as Springs**  
  Covalent bonds between consecutive nucleotides (within a strand) are modeled as linear springs using Hooke’s law. This maintains chain continuity and spacing.

- **Strand Pairing as Springs**  
  Complementary strand connections (e.g., A-T, C-G) are also modeled as springs, enforcing base-pairing constraints and double-helix tendencies.

- **Electrostatic Repulsion**  
  Each nucleotide is treated as having a small partial charge. All nucleotides repel each other based on Coulomb’s law, promoting realistic spacing and discouraging overlaps. The repulsion is computed pairwise across the entire strand.

This hybrid of deterministic and stochastic behavior allows the system to simulate biologically plausible folding patterns with low computational cost. The forces are implemented with tunable constants and are designed to support future upgrades to full 3D force interactions or energy minimization.

---

## 🚀 Features

- Procedural generation of nucleotide positions
- Configurable strand length, noise intensity, and direction
- Modular codebase for easy extension
- Visualization using Matplotlib
