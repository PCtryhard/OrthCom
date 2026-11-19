# OrthCom

Designs a protein pocket for the non-canonical amino acid L-2-naphthylalanine, and
measures whether it binds 2-Nal more tightly than the canonical aromatics 2-Nal
resembles.

## The question

Can a designed pocket bind a synthetic non-canonical amino acid (ncAA) with a
measurable affinity gap over its closest canonical look-alike?

Binding alone settles nothing. A pocket that grips 2-Nal tightly but grips tryptophan
just as tightly has discriminated nothing. The win condition is the gap, and the answer
can be no. For my first target it was no, and I can now say why in kcal/mol.

## Results

**One fluorine is not a discriminable feature.** My first target was
L-4-fluorophenylalanine (L4F) against phenylalanine. Running LigandMPNN on the docked
pocket with ligand atom context on, sequence choice moved between L4F and PHE by a
Jensen-Shannon divergence of 0.0034 bits mean at fluorine-proximal positions, 0.0099
max. Turning whole-ligand context on or off moves the same positions by up to 0.60
bits, roughly 180 times more. At the ring-adjacent position the model picked Phe
deterministically with the fluorine and without it.

**Why, in volume.** Each anchor superposed on the canonical look-alike it has to be
told apart from, van der Waals volume on a 0.15 A grid:

| anchor | vdW volume | surplus over PHE | PHE-shaped void | reach from CB |
|---|---|---|---|---|
| PHE | 97.6 A^3 | | | 4.28 A |
| L4F | 104.3 A^3 | 6.7 A^3 | 0.0 A^3 | 5.62 A |
| 2-Nal | 136.6 A^3 | 39.9 A^3 | 1.0 A^3 | 6.46 A |

The zero in the void column is the whole explanation. PHE fits entirely inside the L4F
envelope, so a pocket complementary to L4F admits PHE with no steric penalty at all.
The only signal available was a favourable contact with the fluorine itself, worth
about 0.2 to 0.3 kcal/mol once 6.7 A^3 is converted to buried nonpolar surface, or
about 1.5x in Kd. Nothing in the pipeline resolves that, so LigandMPNN was not
underperforming. The signal was not there.

2-Nal carries 6.0x the surplus volume and reaches 2.2 A further from CB, order 1.2
kcal/mol, or about 10x in Kd. That is why the target changed.
`python src/analyze_shape_gap.py` reproduces the table from the committed anchors.

**Against tryptophan there is finally something to grip.** Trp is the only
fused-bicyclic canonical aromatic, so it is the real cross-reactivity threat rather than
PHE. Comparing it fairly means giving it its best shot: I fix the attachment at CA, CB
and the ipso carbon, then search the chi2 dihedral for the Trp rotamer that maximises
overlap with 2-Nal. At that best-case rotamer, chi2 180 degrees:

| quantity | value | what it means |
|---|---|---|
| shared volume | 103.9 A^3, 76% of 2-Nal | Trp is a good mimic, which is why this test matters |
| Trp surplus | 19.7 A^3 | volume a 2-Nal pocket does not provide, so Trp has to clash or distort |
| 2-Nal void | 32.5 A^3 | contact Trp cannot make |
| indole NE1 | 0.56 A from the nearest 2-Nal atom, inside the envelope | the N-H donor lands where the pocket will have put apolar surface |

Three independent levers, and all three favour 2-Nal. The 19.7 A^3 surplus is the first
genuine steric term in this project, against 0.0 A^3 for PHE. The 32.5 A^3 void is lost
contact worth roughly 0.8 kcal/mol. And because NE1 sits 0.56 A from a naphthalene
carbon position, a pocket complementary to 2-Nal necessarily buries Trp's N-H against
apolar surface unsatisfied, which costs 1 to 3 kcal/mol and does not depend on volume at
all. That last term is the one the pocket design objective should be written around.
Across the chi2 scan the overlap ranges 72.6 to 103.9 A^3, so I quote the best case
throughout.

**No pocket exists yet.** RFdiffusionAA has not been run, so there is no designed
backbone and no affinity gap measured for either target. What is validated is the
interface instrument described below.

## Method

The design is conditioned on the actual ncAA at every learnable step, so the models see
the atoms that are supposed to be discriminated.

| stage | tool | scripts | state |
|---|---|---|---|
| ncAA anchor | RCSB CCD ideal coords | `src/prep_ligand.py`, `src/prep_ligand_2nal.py`, `src/prep_ligand_trp.py` | L4F, PHE, 2-Nal, Trp built and hash-pinned |
| shape analysis | grid vdW volume | `src/analyze_shape_gap.py` | done, CPU |
| pocket backbone | RFdiffusionAA | `src/prep_rfdiffusion_aa_input.py`, `cloud/rfdiffusion_aa/` | 2-Nal input staged, cloud only, not run |
| design filter | N1 burial and acceptor check | `src/check_n1_pocket.py` | ready, needs a design to run on |
| pocket sequence | LigandMPNN | `src/prep_ligandmpnn_input.py`, `src/analyze_ligandmpnn_dryrun.py` | installed, dry-run done |
| interface validation | ColabFold AF2-multimer | `src/prep_multimer.py`, `src/parse_interface.py` | working |
| affinity gap | Psi4, RESP, AmberTools, MD | `src/generate_conformer.py`, `src/run_psi4_dft.py` | partial, DFT SCF only |

Anchors use a side-chain representation: side chain plus the CA cap, no free-terminus
oxygens, no hydrogens. In the real peptide the ncAA sits mid-chain, so it has no free
carboxylate, and a buried charged carboxylate would dominate the energetics over the
shape difference being tested. `src/prep_ligand_2nal.py` builds 2-Nal from CCD entry
NAL and checks 12 carbons, 0 oxygens, 11 ring bonds and fusion carbons at C4A and C8A
before writing.

Interfaces are scored by co-folding the designed pocket with the frozen peptide as a
two-chain complex and reading cross-chain metrics, not monomer fold quality.
`parse_interface.py` ranks designs by interface PAE, the mean of the two cross-chain
PAE blocks, with iPTM as tie-break. On an 8 GB GPU this only runs stably under the JAX
`platform` allocator; the BFC allocator fragments and segfaults at recycle 2.

## Limitations

The gap against PHE is lost contact, not steric exclusion. The void column shows PHE
fits inside the 2-Nal envelope, so a 2-Nal pocket does not exclude PHE, it just buries
less of it. Only the Trp comparison has a steric term, and even there the 19.7 A^3
surplus assumes a pocket moulded exactly to the 2-Nal envelope with no slack. Real
pockets breathe and Trp can pick a different rotamer, so treat that number as an upper
bound on the steric penalty rather than a prediction.

What has been measured against Trp is shape, not affinity. The three levers above are
geometric arguments about what a pocket could read, not a computed binding free energy,
and no pocket exists yet to test them on. The 2-Nal-vs-Trp affinity gap remains the
single result most likely to falsify this project.

The volume numbers are static van der Waals volumes. They carry no entropy, no
desolvation and no induced fit, so treat the kcal/mol conversions as order-of-magnitude
only. The buried-donor estimate of 1 to 3 kcal/mol is a literature range for
unsatisfied buried donors generally, not a calculation on this system. 2-Nal is held at
its CCD ideal chi2 throughout, since that is the conformer the anchor file carries, so
the comparison does not explore 2-Nal's own rotamer freedom. The affinity gap itself is
unmeasured: the physics path runs DFT single points but does not export the ESP grid or
run the AmberTools RESP step.

Designed proteins have no MSA, so AF2 runs with `--msa-mode single_sequence` and its
confidence should be read accordingly. Everything is in silico. One scaffold (PDB 4QAE)
supplies the carrier protein and the frozen peptide chain, so none of this is varied.

## How to run

```bash
pip install -r requirements.txt

python src/prep_ligand.py              # -> data/raw/ligand_L4F.pdb, ligand_PHE.pdb
python src/prep_ligand_2nal.py         # -> data/raw/ligand_2NAL.pdb, verifies sha1
python src/prep_ligand_trp.py          # -> data/raw/ligand_TRP.pdb, verifies sha1
python src/analyze_shape_gap.py        # both tables above, CPU, no weights, ~1 min
```

The pocket design step needs a GPU with at least 16 GB and Apptainer, so it is cloud
only. `cloud/rfdiffusion_aa/README.md` is a self-contained package with the setup and the
exact run command, and it now carries the 2-Nal input alongside the L4F baseline:

```bash
python src/prep_rfdiffusion_aa_input.py              # -> cloud/.../input_2NAL_pocket.pdb
python src/prep_rfdiffusion_aa_input.py --residue L4F  # regenerates the baseline input
```

After a backbone exists, the N1 criterion is applied as a filter rather than as a design
objective, because RFdiffusionAA conditions on ligand atoms and has no term for it:

```bash
python src/check_n1_pocket.py --pocket output/2nal_pocket/des_0.pdb
```

It requires NE1 to be buried with no hydrogen-bond acceptor within 3.5 A, and exits
non-zero on failure so it can gate a loop over many designs. Then LigandMPNN for
sequence, the AF2 co-fold instrument, and the physics affinity gap against Trp and PHE.

## Notes

`PROGRESS_LOG.md` is my working log: what I tried, what broke, and why I changed
direction. I wrote it as I went, so the early entries describe a target I later dropped,
and I have left them there rather than tidying them away.

Related: [orthcom-immuno](https://github.com/PCtryhard/orthcom-immuno) applies the same
verified ncAA anchors and physics path to a different question, whether non-canonical
residues change MHC-I antigen presentation, and has measured results.
