#!/usr/bin/env python3
"""
Tests verifying that the INDELIBLE workflow correctly handles
CRAM + CRAM.CRAI inputs for individual samples, mother-only,
father-only, and trio family structures.

These tests exercise:
  1. CRAM filename parsing and sample-ID extraction (mirroring the
     Nextflow fromFilePairs logic in main.nf).
  2. Correct grouping into crams / cram_mom / cram_dad / cram_trios
     channels based on which family members are present.
  3. End-to-end TSV → filter → VCF conversion for each family type,
     including chromosome-X entries from the actual INDELIBLE output.

Nextflow channel logic being replicated (from main.nf):

  // ch_crams – proband CRAMs only
  Channel.fromFilePairs(['*{.cram,.cram.crai}'])
      .map { it -> [ it[0][0..-6], it[1][0], it[1][1] ] }   // strip _01_1
      .filter { it -> it[1] =~ '_01_1' }

  // ch_cram_trios – complete trios (6 files)
  Channel.fromFilePairs(['*{_01_1,_02_2,_03_3}*'], size: 6)

  // ch_cram_mom – child + mother only (4 files, last file is _02_2.cram*)
  Channel.fromFilePairs(['*{_01_1,_02_2,_03_3}*'], size: -1)
      .filter { it -> it[1].size() == 4 && it[1][3] =~ '02_2.cram' }

  // ch_cram_dad – child + father only (4 files, last file is _03_3.cram*)
  Channel.fromFilePairs(['*{_01_1,_02_2,_03_3}*'], size: -1)
      .filter { it -> it[1].size() == 4 && it[1][3] =~ '03_3.cram' }
"""

import os
import sys
from collections import defaultdict

import pytest

BIN_DIR = os.path.join(os.path.dirname(__file__), '..', 'bin')
sys.path.insert(0, os.path.abspath(BIN_DIR))

# ---------------------------------------------------------------------------
# Channel construction helpers (Python simulation of Nextflow fromFilePairs)
# ---------------------------------------------------------------------------

def _sample_id_from_stem(stem):
    """
    Groovy:  it[0][0..-6]
    Strips the last 5 characters from the fromFilePairs key, which removes
    the role suffix (_01_1 / _02_2 / _03_3, each exactly 5 chars).
    """
    return stem[:-5]


def _group_cram_files(cram_dir):
    """
    Simulate Nextflow's INDELIBLE channel construction from main.nf.

    Given a directory containing .cram and .cram.crai files following the
    naming convention  <family_id>_{01_1|02_2|03_3}.cram[.crai], return:

      crams      – list of (family_id, cram, crai)         [proband only]
      cram_trios – list of (family_id, c_cram, c_crai, m_cram, m_crai, d_cram, d_crai)
      cram_mom   – list of (family_id, c_cram, c_crai, m_cram, m_crai)
      cram_dad   – list of (family_id, c_cram, c_crai, d_cram, d_crai)
    """
    # --- Step 1: build CRAM/CRAI pairs (mirrors fromFilePairs *{.cram,.cram.crai}) ---
    pairs = {}   # stem -> [cram_path, crai_path]
    for fname in sorted(os.listdir(cram_dir)):
        fpath = os.path.join(cram_dir, fname)
        if fname.endswith('.cram.crai'):
            stem = fname[:-len('.cram.crai')]
            pairs.setdefault(stem, [None, None])[1] = fpath
        elif fname.endswith('.cram'):
            stem = fname[:-len('.cram')]
            pairs.setdefault(stem, [None, None])[0] = fpath

    # keep only complete pairs
    complete = {s: v for s, v in pairs.items() if v[0] and v[1]}

    # --- Step 2: ch_crams – proband (filter to _01_1 stems, strip suffix) ---
    crams = []
    for stem, (cram, crai) in sorted(complete.items()):
        if '_01_1' in stem:
            family_id = _sample_id_from_stem(stem)
            crams.append((family_id, cram, crai))

    # --- Step 3: group by family for trio/mom/dad channels ---
    # fromFilePairs with *{_01_1,_02_2,_03_3}* groups by the common prefix
    family_files = defaultdict(list)
    for stem, (cram, crai) in sorted(complete.items()):
        for role in ('_01_1', '_02_2', '_03_3'):
            if role in stem:
                prefix = stem[:stem.index(role)]
                family_files[prefix].extend([cram, crai])
                break

    cram_trios, cram_mom, cram_dad = [], [], []
    for family_id, files in sorted(family_files.items()):
        files = sorted(files)       # lexicographic – mirrors Nextflow
        if len(files) == 6:
            # size: 6 → complete trio: [c.cram, c.crai, m.cram, m.crai, d.cram, d.crai]
            cram_trios.append((family_id, *files))
        elif len(files) == 4:
            # size: -1 filtered to 4 – mom or dad
            # it[1][3] =~ '02_2.cram'  →  fourth file contains '02_2.cram'
            if '02_2.cram' in files[3]:
                cram_mom.append((family_id, *files))
            # it[1][3] =~ '03_3.cram'  →  fourth file contains '03_3.cram'
            elif '03_3.cram' in files[3]:
                cram_dad.append((family_id, *files))

    return crams, cram_trios, cram_mom, cram_dad


# ---------------------------------------------------------------------------
# Fixtures: minimal fake CRAM / CRAI files
# ---------------------------------------------------------------------------

def _touch(path):
    """Create an empty file."""
    open(path, 'w').close()


@pytest.fixture()
def cram_dir_solo(tmp_path):
    """A directory with a single proband CRAM+CRAI (no parents)."""
    d = tmp_path / "solo"
    d.mkdir()
    _touch(d / "FAMILY1_01_1.cram")
    _touch(d / "FAMILY1_01_1.cram.crai")
    return str(d)


@pytest.fixture()
def cram_dir_mom(tmp_path):
    """A directory with proband + mother CRAMs."""
    d = tmp_path / "mom"
    d.mkdir()
    for f in ["FAMILY2_01_1.cram", "FAMILY2_01_1.cram.crai",
              "FAMILY2_02_2.cram", "FAMILY2_02_2.cram.crai"]:
        _touch(d / f)
    return str(d)


@pytest.fixture()
def cram_dir_dad(tmp_path):
    """A directory with proband + father CRAMs."""
    d = tmp_path / "dad"
    d.mkdir()
    for f in ["FAMILY3_01_1.cram", "FAMILY3_01_1.cram.crai",
              "FAMILY3_03_3.cram", "FAMILY3_03_3.cram.crai"]:
        _touch(d / f)
    return str(d)


@pytest.fixture()
def cram_dir_trio(tmp_path):
    """A directory with a complete trio (proband + mother + father)."""
    d = tmp_path / "trio"
    d.mkdir()
    for f in ["FAMILY4_01_1.cram", "FAMILY4_01_1.cram.crai",
              "FAMILY4_02_2.cram", "FAMILY4_02_2.cram.crai",
              "FAMILY4_03_3.cram", "FAMILY4_03_3.cram.crai"]:
        _touch(d / f)
    return str(d)


@pytest.fixture()
def cram_dir_multi(tmp_path):
    """A directory with mixed family structures: solo, mom, dad, trio."""
    d = tmp_path / "multi"
    d.mkdir()
    # Solo proband
    for f in ["SOLO_01_1.cram", "SOLO_01_1.cram.crai"]:
        _touch(d / f)
    # Mom family
    for f in ["MOM_01_1.cram", "MOM_01_1.cram.crai",
              "MOM_02_2.cram", "MOM_02_2.cram.crai"]:
        _touch(d / f)
    # Dad family
    for f in ["DAD_01_1.cram", "DAD_01_1.cram.crai",
              "DAD_03_3.cram", "DAD_03_3.cram.crai"]:
        _touch(d / f)
    # Trio family
    for f in ["TRIO_01_1.cram", "TRIO_01_1.cram.crai",
              "TRIO_02_2.cram", "TRIO_02_2.cram.crai",
              "TRIO_03_3.cram", "TRIO_03_3.cram.crai"]:
        _touch(d / f)
    return str(d)


# ---------------------------------------------------------------------------
# TSV fixtures using actual INDELIBLE output format (from problem statement)
# ---------------------------------------------------------------------------

# The 37-column header matching the real INDELIBLE annotated TSV output
_INDELIBLE_HEADER = (
    "chrom\tposition\tcoverage\tinsertion_context\tdeletion_context\t"
    "sr_total\tsr_total_long\tsr_total_short\tsr_long_5\tsr_short_5\t"
    "sr_long_3\tsr_short_3\tsr_entropy\tcontext_entropy\t"
    "entropy_upstream\tentropy_downstream\tsr_sw_similarity\t"
    "avg_avg_sr_qual\tavg_mapq\tseq_longest\tpredicted\tprob_N\tprob_Y\t"
    "ddg2p\thgnc\thgnc_constrained\texonic\ttranscripts\texon_numbers\t"
    "maf\tblast_hit\tblast_strand\tblast_identity\tblast_dist\t"
    "blast_hgnc\tblast_hgnc_constrained\tblast_ddg2p"
)

# Two verbatim rows from the problem statement
_INDELIBLE_ROW_N = (
    "X\t153296122\t19\t0\t0\t3\t2\t1\t2\t1\t0\t0\t"
    "1.71473664918\t1.65515020808\t1.70600757931\t1.59546184424\t"
    "1.0\t32.835\t29.0\t"
    "ATTTCCCTTTTTACTTAAGATCATTTCAGT\t"
    "N\t0.911525\t0.088475\t"
    "MECP2\tMECP2\t\tFalse\t\t\t\t"
    "2:1610630-1610601\tminus\t100.0\tother_chrom\tNA\tNA\tNA"
)

_INDELIBLE_ROW_Y = (
    "X\t153296083\t31\t0\t0\t20\t17\t3\t0\t0\t17\t3\t"
    "1.99175543619\t1.64419231285\t1.44115352788\t1.73935387217\t"
    "1.0\t32.5648148148\t29.0\t"
    "TCATATTGACCTTGGCTGCAGTGCCGGACAAGGAGCCTGTAAGGTGCAGTCACTC\t"
    "Y\t0.23161\t0.76839\t"
    "MECP2\tMECP2\t\tFalse\t\t\t\t"
    "2:1679342-1679288\tminus\t100.0\tother_chrom\tNA\tNA\tNA"
)


@pytest.fixture()
def indelible_chrX_tsv(tmp_path):
    """Annotated TSV with the two real chromosome-X rows from the problem statement."""
    tsv = tmp_path / "MECP2_SAMPLE.counts.scored.annotated"
    tsv.write_text(
        _INDELIBLE_HEADER + "\n" +
        _INDELIBLE_ROW_N + "\n" +
        _INDELIBLE_ROW_Y + "\n"
    )
    return str(tsv)


@pytest.fixture()
def fai_with_chrX(tmp_path):
    """Minimal FAI file that includes chrX."""
    fai = tmp_path / "ref.fa.fai"
    fai.write_text(
        "chr1\t248956422\t112\t70\t71\n"
        "chrX\t156040895\t2702094775\t70\t71\n"
        "chrY\t57227415\t2860899383\t70\t71\n"
    )
    return str(fai)


# ---------------------------------------------------------------------------
# Tests: sample-ID extraction from CRAM filenames
# ---------------------------------------------------------------------------

class TestCramSampleIdExtraction:
    """Sample ID extraction mirrors Nextflow's it[0][0..-6] on the fromFilePairs key."""

    def test_proband_suffix_stripped(self):
        """_01_1 (5 chars) is stripped by [0..-6] → base family ID."""
        assert _sample_id_from_stem("FAMILY_01_1") == "FAMILY"

    def test_mother_suffix_stripped(self):
        """_02_2 (5 chars) is stripped by [0..-6]."""
        assert _sample_id_from_stem("FAMILY_02_2") == "FAMILY"

    def test_father_suffix_stripped(self):
        """_03_3 (5 chars) is stripped by [0..-6]."""
        assert _sample_id_from_stem("FAMILY_03_3") == "FAMILY"

    def test_complex_family_id(self):
        """Multi-part family IDs are preserved; only the trailing role suffix is removed."""
        assert _sample_id_from_stem("FAM123_PROBAND_01_1") == "FAM123_PROBAND"

    def test_strip_length_is_five(self):
        """The trailing suffix removed is exactly 5 characters (_01_1 etc.)."""
        stem = "X_01_1"
        assert len(stem) - len(_sample_id_from_stem(stem)) == 5


# ---------------------------------------------------------------------------
# Tests: channel construction – proband-only (ch_crams)
# ---------------------------------------------------------------------------

class TestCramChannelSolo:
    """ch_crams should contain only the proband CRAM+CRAI."""

    def test_proband_cram_present(self, cram_dir_solo):
        crams, _, _, _ = _group_cram_files(cram_dir_solo)
        assert len(crams) == 1

    def test_sample_id_is_family_prefix(self, cram_dir_solo):
        crams, _, _, _ = _group_cram_files(cram_dir_solo)
        family_id, cram, crai = crams[0]
        assert family_id == "FAMILY1"

    def test_cram_file_has_cram_extension(self, cram_dir_solo):
        crams, _, _, _ = _group_cram_files(cram_dir_solo)
        _, cram, _ = crams[0]
        assert cram.endswith('.cram') and not cram.endswith('.cram.crai')

    def test_crai_file_has_crai_extension(self, cram_dir_solo):
        crams, _, _, _ = _group_cram_files(cram_dir_solo)
        _, _, crai = crams[0]
        assert crai.endswith('.cram.crai')

    def test_no_trio_mom_dad_channels(self, cram_dir_solo):
        _, trios, moms, dads = _group_cram_files(cram_dir_solo)
        assert trios == []
        assert moms == []
        assert dads == []


# ---------------------------------------------------------------------------
# Tests: channel construction – proband + mother (ch_cram_mom)
# ---------------------------------------------------------------------------

class TestCramChannelMom:
    """ch_cram_mom should be populated; ch_cram_trios and ch_cram_dad empty."""

    def test_mom_channel_populated(self, cram_dir_mom):
        _, _, moms, _ = _group_cram_files(cram_dir_mom)
        assert len(moms) == 1

    def test_mom_channel_has_five_elements(self, cram_dir_mom):
        """Tuple: (family_id, child_cram, child_crai, mom_cram, mom_crai)."""
        _, _, moms, _ = _group_cram_files(cram_dir_mom)
        assert len(moms[0]) == 5

    def test_mom_child_cram_contains_01_1(self, cram_dir_mom):
        _, _, moms, _ = _group_cram_files(cram_dir_mom)
        _, child_cram, _, _, _ = moms[0]
        assert '_01_1' in child_cram

    def test_mom_cram_contains_02_2(self, cram_dir_mom):
        _, _, moms, _ = _group_cram_files(cram_dir_mom)
        _, _, _, mom_cram, _ = moms[0]
        assert '_02_2' in mom_cram

    def test_crams_channel_still_has_proband(self, cram_dir_mom):
        crams, _, _, _ = _group_cram_files(cram_dir_mom)
        assert len(crams) == 1

    def test_trio_and_dad_channels_empty(self, cram_dir_mom):
        _, trios, _, dads = _group_cram_files(cram_dir_mom)
        assert trios == []
        assert dads == []


# ---------------------------------------------------------------------------
# Tests: channel construction – proband + father (ch_cram_dad)
# ---------------------------------------------------------------------------

class TestCramChannelDad:
    """ch_cram_dad should be populated; ch_cram_trios and ch_cram_mom empty."""

    def test_dad_channel_populated(self, cram_dir_dad):
        _, _, _, dads = _group_cram_files(cram_dir_dad)
        assert len(dads) == 1

    def test_dad_channel_has_five_elements(self, cram_dir_dad):
        """Tuple: (family_id, child_cram, child_crai, dad_cram, dad_crai)."""
        _, _, _, dads = _group_cram_files(cram_dir_dad)
        assert len(dads[0]) == 5

    def test_dad_child_cram_contains_01_1(self, cram_dir_dad):
        _, _, _, dads = _group_cram_files(cram_dir_dad)
        _, child_cram, _, _, _ = dads[0]
        assert '_01_1' in child_cram

    def test_dad_cram_contains_03_3(self, cram_dir_dad):
        _, _, _, dads = _group_cram_files(cram_dir_dad)
        _, _, _, dad_cram, _ = dads[0]
        assert '_03_3' in dad_cram

    def test_crams_channel_still_has_proband(self, cram_dir_dad):
        crams, _, _, _ = _group_cram_files(cram_dir_dad)
        assert len(crams) == 1

    def test_trio_and_mom_channels_empty(self, cram_dir_dad):
        _, trios, moms, _ = _group_cram_files(cram_dir_dad)
        assert trios == []
        assert moms == []


# ---------------------------------------------------------------------------
# Tests: channel construction – complete trio (ch_cram_trios)
# ---------------------------------------------------------------------------

class TestCramChannelTrio:
    """ch_cram_trios should be populated; ch_cram_mom and ch_cram_dad empty."""

    def test_trio_channel_populated(self, cram_dir_trio):
        _, trios, _, _ = _group_cram_files(cram_dir_trio)
        assert len(trios) == 1

    def test_trio_channel_has_seven_elements(self, cram_dir_trio):
        """Tuple: (family_id, c_cram, c_crai, m_cram, m_crai, d_cram, d_crai)."""
        _, trios, _, _ = _group_cram_files(cram_dir_trio)
        assert len(trios[0]) == 7

    def test_trio_file_order_child_mom_dad(self, cram_dir_trio):
        """Lexicographic sort places _01_1 < _02_2 < _03_3."""
        _, trios, _, _ = _group_cram_files(cram_dir_trio)
        _, c_cram, c_crai, m_cram, m_crai, d_cram, d_crai = trios[0]
        assert '_01_1' in c_cram
        assert '_01_1' in c_crai
        assert '_02_2' in m_cram
        assert '_02_2' in m_crai
        assert '_03_3' in d_cram
        assert '_03_3' in d_crai

    def test_crams_channel_still_has_proband(self, cram_dir_trio):
        crams, _, _, _ = _group_cram_files(cram_dir_trio)
        assert len(crams) == 1

    def test_mom_and_dad_channels_empty(self, cram_dir_trio):
        _, _, moms, dads = _group_cram_files(cram_dir_trio)
        assert moms == []
        assert dads == []


# ---------------------------------------------------------------------------
# Tests: multi-family directory (all structures present simultaneously)
# ---------------------------------------------------------------------------

class TestCramChannelMultiFamily:
    """Multiple families with different structures co-exist in one directory."""

    def test_four_probands_in_crams_channel(self, cram_dir_multi):
        crams, _, _, _ = _group_cram_files(cram_dir_multi)
        assert len(crams) == 4

    def test_one_trio(self, cram_dir_multi):
        _, trios, _, _ = _group_cram_files(cram_dir_multi)
        assert len(trios) == 1

    def test_one_mom_family(self, cram_dir_multi):
        _, _, moms, _ = _group_cram_files(cram_dir_multi)
        assert len(moms) == 1

    def test_one_dad_family(self, cram_dir_multi):
        _, _, _, dads = _group_cram_files(cram_dir_multi)
        assert len(dads) == 1

    def test_all_proband_sample_ids_are_correct(self, cram_dir_multi):
        crams, _, _, _ = _group_cram_files(cram_dir_multi)
        ids = {t[0] for t in crams}
        assert ids == {"SOLO", "MOM", "DAD", "TRIO"}


# ---------------------------------------------------------------------------
# Tests: VCF conversion with real INDELIBLE TSV (chromosome X data)
# ---------------------------------------------------------------------------

class TestIndelibleChrXConversion:
    """
    VCF conversion of the actual two INDELIBLE rows from the problem statement:
      - Row 1: chrX:153296122  predicted=N  (LowQuality)
      - Row 2: chrX:153296083  predicted=Y  (PASS)
    """

    def test_vcf_created_from_chrX_tsv(self, indelible_chrX_tsv, fai_with_chrX, tmp_path):
        """VCF file is produced without errors from chrX INDELIBLE data."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_chrX_tsv, out, "MECP2_SAMPLE", fai_with_chrX)
        assert os.path.isfile(os.path.join(out, "MECP2_SAMPLE_INDELIBLE_output.vcf"))

    def test_two_records_produced(self, indelible_chrX_tsv, fai_with_chrX, tmp_path):
        """One VCF record is produced for each TSV data row."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_chrX_tsv, out, "MECP2_SAMPLE", fai_with_chrX)

        vcf = os.path.join(out, "MECP2_SAMPLE_INDELIBLE_output.vcf")
        data = [l for l in open(vcf) if not l.startswith('#') and l.strip()]
        assert len(data) == 2, f"Expected 2 records, got {len(data)}"

    def test_chromosome_x_prefixed_as_chrX(self, indelible_chrX_tsv, fai_with_chrX, tmp_path):
        """'X' in the TSV chrom column is normalised to 'chrX' in the VCF."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_chrX_tsv, out, "MECP2_SAMPLE", fai_with_chrX)

        vcf = os.path.join(out, "MECP2_SAMPLE_INDELIBLE_output.vcf")
        data = [l for l in open(vcf) if not l.startswith('#') and l.strip()]
        for record in data:
            record_field = {record.split('\t')[0]}
            assert record.startswith("chrX\t"), (
                f"Expected chrX prefix, got: ${record_field}"
            )

    def test_predicted_N_is_lowquality(self, indelible_chrX_tsv, fai_with_chrX, tmp_path):
        """Row with predicted=N gets FILTER=LowQuality."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_chrX_tsv, out, "MECP2_SAMPLE", fai_with_chrX)

        vcf = os.path.join(out, "MECP2_SAMPLE_INDELIBLE_output.vcf")
        data = [l for l in open(vcf) if not l.startswith('#') and l.strip()]
        # Row 1 (pos 153296122) → predicted=N
        record_n = next(r for r in data if r.split('\t')[1] == '153296122')
        record_n_field = {record_n.split('\t')[6]}
        assert record_n.split('\t')[6] == 'LowQuality', (
            f"Expected LowQuality, got ${record_n_field}"
        )

    def test_predicted_Y_is_pass(self, indelible_chrX_tsv, fai_with_chrX, tmp_path):
        """Row with predicted=Y gets FILTER=PASS."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_chrX_tsv, out, "MECP2_SAMPLE", fai_with_chrX)

        vcf = os.path.join(out, "MECP2_SAMPLE_INDELIBLE_output.vcf")
        data = [l for l in open(vcf) if not l.startswith('#') and l.strip()]
        # Row 2 (pos 153296083) → predicted=Y
        record_y = next(r for r in data if r.split('\t')[1] == '153296083')
        record_y_field = {record_y.split('\t')[6]}
        assert record_y.split('\t')[6] == 'PASS', (
            f"Expected PASS, got ${record_y_field}"
        )

    def test_hgnc_mecp2_in_info(self, indelible_chrX_tsv, fai_with_chrX, tmp_path):
        """HGNC=MECP2 appears in the INFO field of both records."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_chrX_tsv, out, "MECP2_SAMPLE", fai_with_chrX)

        vcf = os.path.join(out, "MECP2_SAMPLE_INDELIBLE_output.vcf")
        data = [l for l in open(vcf) if not l.startswith('#') and l.strip()]
        for record in data:
            info = record.split('\t')[7]
            assert "HGNC=MECP2" in info, f"HGNC=MECP2 not found in INFO: {info}"

    def test_seq_longest_in_info(self, indelible_chrX_tsv, fai_with_chrX, tmp_path):
        """SEQ field is populated from the seq_longest column."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_chrX_tsv, out, "MECP2_SAMPLE", fai_with_chrX)

        vcf = os.path.join(out, "MECP2_SAMPLE_INDELIBLE_output.vcf")
        data = [l for l in open(vcf) if not l.startswith('#') and l.strip()]
        for record in data:
            info = record.split('\t')[7]
            assert "SEQ=" in info, f"SEQ field missing from INFO: {info}"

    def test_exonic_false_in_info(self, indelible_chrX_tsv, fai_with_chrX, tmp_path):
        """EXONIC=False appears in the INFO field (actual INDELIBLE value)."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(indelible_chrX_tsv, out, "MECP2_SAMPLE", fai_with_chrX)

        vcf = os.path.join(out, "MECP2_SAMPLE_INDELIBLE_output.vcf")
        data = [l for l in open(vcf) if not l.startswith('#') and l.strip()]
        for record in data:
            info = record.split('\t')[7]
            assert "EXONIC=False" in info, f"EXONIC=False not found in INFO: {info}"


# ---------------------------------------------------------------------------
# Tests: per-family-type VCF conversion (proband / mother / father / trio)
# ---------------------------------------------------------------------------

class TestIndeliblePerFamilyTypeVcf:
    """
    Each family member type (proband _01_1, mother _02_2, father _03_3)
    produces its own VCF through the standard TSV→filter→VCF pipeline.
    """

    def _make_tsv(self, tmp_path, sample, chrom="chr1", pos="1000",
                  predicted="Y", prob_Y="0.9"):
        """Create a minimal 37-column INDELIBLE annotated TSV for a sample."""
        prob_N = round(1.0 - float(prob_Y), 6)
        row = (
            f"{chrom}\t{pos}\t30\t0\t0\t"
            "15\t10\t5\t7\t3\t3\t2\t"
            "1.8\t1.5\t1.6\t1.4\t1.0\t33.0\t30.0\t"
            f"ACGTACGT\t{predicted}\t{prob_N}\t{prob_Y}\t"
            "NA\tGENE1\tNA\tFalse\tENST001\t2\t0.005\t"
            "NA\tNA\tNA\tNA\tNA\tNA\tNA"
        )
        tsv = tmp_path / f"{sample}.counts.scored.annotated"
        tsv.write_text(_INDELIBLE_HEADER + "\n" + row + "\n")
        return str(tsv)

    def test_proband_vcf_created(self, tmp_path, fai_with_chrX):
        import indelible_tsv_to_vcf as mod
        sample = "FAMILY_01_1"
        tsv = self._make_tsv(tmp_path, sample)
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(tsv, out, sample, fai_with_chrX)
        assert os.path.isfile(os.path.join(out, f"{sample}_INDELIBLE_output.vcf"))

    def test_mother_vcf_created(self, tmp_path, fai_with_chrX):
        import indelible_tsv_to_vcf as mod
        sample = "FAMILY_02_2"
        tsv = self._make_tsv(tmp_path, sample, predicted="N", prob_Y="0.1")
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(tsv, out, sample, fai_with_chrX)
        assert os.path.isfile(os.path.join(out, f"{sample}_INDELIBLE_output.vcf"))

    def test_father_vcf_created(self, tmp_path, fai_with_chrX):
        import indelible_tsv_to_vcf as mod
        sample = "FAMILY_03_3"
        tsv = self._make_tsv(tmp_path, sample, predicted="N", prob_Y="0.2")
        out = str(tmp_path / "out")
        mod.convert_indelible_tsv_to_vcf(tsv, out, sample, fai_with_chrX)
        assert os.path.isfile(os.path.join(out, f"{sample}_INDELIBLE_output.vcf"))

    def test_all_three_vcfs_created_for_trio(self, tmp_path, fai_with_chrX):
        """Running conversion for child + mother + father produces three VCFs."""
        import indelible_tsv_to_vcf as mod
        out = str(tmp_path / "out")
        for sample, pred, prob in [
            ("FAM_01_1", "Y", "0.9"),
            ("FAM_02_2", "N", "0.15"),
            ("FAM_03_3", "N", "0.20"),
        ]:
            tsv = self._make_tsv(tmp_path, sample,
                                 predicted=pred, prob_Y=prob)
            mod.convert_indelible_tsv_to_vcf(tsv, out, sample, fai_with_chrX)

        for sample in ["FAM_01_1", "FAM_02_2", "FAM_03_3"]:
            assert os.path.isfile(os.path.join(out, f"{sample}_INDELIBLE_output.vcf")), (
                f"Missing VCF for {sample}"
            )
