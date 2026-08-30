#!/usr/bin/env python3
"""Build sequence deliverables for OVA P3-13R tetramer/multimer candidates.

The script intentionally uses only the Python standard library. It reads the
source XLSX as OOXML, preserves every original workbook member, and appends two
new worksheets containing designs and QC results.
"""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile
import re
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parent
SOURCE_XLSX = ROOT / "OVA WT及OVA-P3-1P突变体序列.xlsx"
OUTPUT_XLSX = ROOT / "OVA_P3-13R_四聚体候选序列.xlsx"
OUTPUT_AA = ROOT / "OVA_P3-13R_四聚体候选_AA.fasta"
OUTPUT_CDS = ROOT / "OVA_P3-13R_四聚体候选_CDS.fasta"
TOP_OUTPUT_AA = ROOT / "OVA-BODY-C4-C3HOT-AI01_AA.fasta"
TOP_OUTPUT_CDS = ROOT / "OVA-BODY-C4-C3HOT-AI01_CDS.fasta"
SS3_OUTPUT_AA = ROOT / "OVA-BODY-C4-SS3-AI03_AA.fasta"
SS3_OUTPUT_CDS = ROOT / "OVA-BODY-C4-SS3-AI03_CDS.fasta"

MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
CT_NS = "http://schemas.openxmlformats.org/package/2006/content-types"
APP_NS = "http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"
VT_NS = "http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes"

ET.register_namespace("", MAIN_NS)
ET.register_namespace("r", REL_NS)


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

RESIDUE_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}

# Deterministic mammalian-biased codons for newly designed peptide modules.
# The source OVA CDS is never reverse-translated: it is copied byte-for-byte
# from the supplied workbook.
PREFERRED_CODON = {
    "A": "GCC", "R": "CGC", "N": "AAC", "D": "GAC", "C": "TGC",
    "E": "GAG", "Q": "CAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
}


def qname(ns: str, local: str) -> str:
    return f"{{{ns}}}{local}"


def translate(dna: str, terminal_stop: bool = False) -> str:
    dna = re.sub(r"\s+", "", dna).upper()
    assert len(dna) % 3 == 0, f"CDS length is not divisible by 3: {len(dna)}"
    aa = "".join(CODON_TABLE[dna[i : i + 3]] for i in range(0, len(dna), 3))
    if terminal_stop:
        assert aa.endswith("*"), "Expected a terminal stop codon"
        assert "*" not in aa[:-1], "Internal stop codon found"
        return aa[:-1]
    assert "*" not in aa, "Unexpected stop codon found"
    return aa


def reverse_translate(aa: str) -> str:
    dna = "".join(PREFERRED_CODON[residue] for residue in aa)
    assert translate(dna) == aa
    return dna


def mutate_source_cds(source_dna: str, changes: dict[int, str]) -> str:
    """Change selected codons while preserving every other source nucleotide."""
    codons = [source_dna[index:index + 3] for index in range(0, len(source_dna), 3)]
    source_aa = translate(source_dna)
    for position, amino_acid in changes.items():
        assert 1 <= position <= len(codons)
        assert source_aa[position - 1] == CODON_TABLE[codons[position - 1]]
        codons[position - 1] = PREFERRED_CODON[amino_acid]
    output = "".join(codons)
    expected = list(source_aa)
    for position, amino_acid in changes.items(): expected[position - 1] = amino_acid
    assert translate(output) == "".join(expected)
    return output


def molecular_weight_kda(aa: str) -> float:
    return (sum(RESIDUE_MASS[x] for x in aa) + 18.0153) / 1000.0


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def read_source_sequences() -> dict[str, str]:
    ns = {"m": MAIN_NS, "r": REL_NS, "pr": PKG_REL_NS}
    with ZipFile(SOURCE_XLSX) as zf:
        shared_root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
        shared = [
            "".join(node.text or "" for node in item.iter(qname(MAIN_NS, "t")))
            for item in shared_root.findall("m:si", ns)
        ]
        sheet_root = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
        cells: dict[str, str] = {}
        for cell in sheet_root.findall(".//m:sheetData/m:row/m:c", ns):
            ref = cell.attrib["r"]
            value = cell.find("m:v", ns)
            if value is None:
                continue
            cells[ref] = shared[int(value.text)] if cell.attrib.get("t") == "s" else value.text or ""
    return {"OVA WT": cells["B1"], "OVA P3-13R": cells["B2"]}


def make_constructs() -> tuple[list[dict[str, object]], dict[str, str], list[tuple[str, str]]]:
    source = read_source_sequences()
    wt_dna = source["OVA WT"]
    p3_dna = source["OVA P3-13R"]
    wt_aa = translate(wt_dna)
    p3_aa = translate(p3_dna)

    assert len(wt_dna) == len(p3_dna) == 1158
    assert len(wt_aa) == len(p3_aa) == 386
    expected_mutations = [(14, "D", "R"), (69, "S", "A"), (77, "S", "A"),
                          (79, "N", "A"), (191, "D", "A"), (192, "E", "F"),
                          (193, "D", "A")]
    observed = [(i + 1, a, b) for i, (a, b) in enumerate(zip(wt_aa, p3_aa)) if a != b]
    assert observed == expected_mutations, f"Unexpected source mutations: {observed}"

    # Exact five-residue linker used in the AF3/ProteinMPNN candidates.
    linker5_aa = "GGGGS"
    linker5_dna = "GGTGGCGGAGGGTCC"
    assert translate(linker5_dna) == linker5_aa

    # GCN4-pLI, PDB 1GCL: RMKQIEDKLEEILSKLYHIENELARIKKLLGER.
    # The DNA below is a generic mammalian-biased reverse translation; OVA itself
    # is preserved exactly as supplied in the source workbook.
    pli_dna = (
        "CGCATGAAGCAGATCGAGGACAAGCTGGAGGAGATCCTGAGCAAGCTG"
        "TACCACATCGAGAACGAGCTGGCCCGCATCAAGAAGCTGCTGGGCGAGCGC"
    )
    pli_aa = translate(pli_dna)
    assert pli_aa == "RMKQIEDKLEEILSKLYHIENELARIKKLLGER"

    ai_c4_04_aa = "NTNVILNRLTTIEAQLVRIELTLARIERLLGVK"
    ai_c4_06_aa = "PTNVITARLTTIEYKLDVIELKLDRIEKLLGVK"
    ai_c4_04_dna = reverse_translate(ai_c4_04_aa)
    ai_c4_06_dna = reverse_translate(ai_c4_06_aa)
    bindcraft_aa = "DAEIQKKIEEIRKKYEKKREEELKKIKEEFQEKEKKGKPTNADYIDYFDKRWEVYRKTNRAYDEEAAKVL"
    bindcraft_dna = reverse_translate(bindcraft_aa)
    body_ss3_changes = {
        93: "E", 94: "C", 96: "C", 97: "N", 99: "Q", 101: "L",
        152: "A", 155: "C", 190: "C", 195: "R", 196: "C", 213: "C",
        237: "D", 239: "Q", 340: "A",
    }
    body_ss3_aa_chars = list(p3_aa)
    for position, amino_acid in body_ss3_changes.items():
        body_ss3_aa_chars[position - 1] = amino_acid
    body_ss3_aa = "".join(body_ss3_aa_chars)
    body_ss3_dna = mutate_source_cds(p3_dna, body_ss3_changes)
    assert translate(body_ss3_dna) == body_ss3_aa
    body_c3hot_changes = {
        90: "E", 93: "E", 94: "C", 96: "C", 97: "N", 99: "Q", 101: "L",
        152: "A", 155: "C", 182: "R", 190: "C", 194: "A", 195: "R",
        196: "C", 210: "L", 213: "C", 214: "T", 237: "D", 340: "A",
    }
    body_c3hot_aa_chars = list(p3_aa)
    for position, amino_acid in body_c3hot_changes.items():
        body_c3hot_aa_chars[position - 1] = amino_acid
    body_c3hot_aa = "".join(body_c3hot_aa_chars)
    body_c3hot_dna = mutate_source_cds(p3_dna, body_c3hot_changes)
    assert translate(body_c3hot_dna) == body_c3hot_aa
    body_p208f_changes = {
        93: "E", 96: "C", 97: "N", 99: "Q", 101: "L", 152: "A",
        153: "A", 155: "C", 190: "C", 195: "R", 196: "C", 208: "F",
        237: "D", 239: "Q", 340: "A",
    }
    body_p208f_aa_chars = list(p3_aa)
    for position, amino_acid in body_p208f_changes.items():
        body_p208f_aa_chars[position - 1] = amino_acid
    body_p208f_aa = "".join(body_p208f_aa_chars)
    body_p208f_dna = mutate_source_cds(p3_dna, body_p208f_changes)
    assert translate(body_p208f_dna) == body_p208f_aa

    def item(
        construct_id: str,
        priority: str,
        purpose: str,
        architecture: str,
        assembly: str,
        aa: str,
        cds: str,
        assembly_factor: int,
        risk: str,
    ) -> dict[str, object]:
        assert cds.endswith("TAA")
        assert translate(cds, terminal_stop=True) == aa
        chain_mw = molecular_weight_kda(aa)
        return {
            "id": construct_id,
            "priority": priority,
            "purpose": purpose,
            "architecture": architecture,
            "assembly": assembly,
            "aa": aa,
            "cds": cds,
            "assembly_factor": assembly_factor,
            "chain_mw": chain_mw,
            "assembly_mw": chain_mw * assembly_factor,
            "risk": risk,
        }

    constructs = [
        item(
            "D0-P3-13R",
            "必做对照",
            "复现实验起点",
            "OVA P3-13R",
            "AF3 未支持同源二聚；按单链/弱缔合对照解释",
            p3_aa,
            p3_dna + "TAA",
            1,
            "25 个 AF3 二聚模型：ipTM 中位 0.17、界面 PAE 22.61 Å；不能称为 AF3 已验证二聚体。",
        ),
        item(
            "D0-WT",
            "必做对照",
            "单体/背景对照",
            "OVA WT",
            "主要按单体解释",
            wt_aa,
            wt_dna + "TAA",
            1,
            "用于判断新增模块与 P3 突变各自的贡献。",
        ),
        item(
            "OVA-BODY-C4-C3HOT-AI01",
            "当前纯本体 C4-forming 首选（链数不特异）",
            "C4/C3 AF3 接触差热点 + ProteinMPNN 多状态负设计 + 三组定向二硫键",
            "386 aa OVA P3-13R 本体点突变；无 linker、binder 或外接四聚螺旋",
            "AF3 稳健四链同源装配；但 C2/C3 也可形成，不能称四聚体特异",
            body_c3hot_aa,
            body_c3hot_dna + "TAA",
            4,
            "八 seed 40/40 高置信四链连通、40/40 无 clash；ipTM/最弱链 0.79/0.78、界面 PAE 1.723 Å；39/40 化学干净完整 12 键网络、0 非预期 S–S。C2/C3/C5/C6 为 0.74/0.75/0.36/0.34；C4 最高且 C5/C6 失败，但 C2/C3 仍强，故只称 C4-forming 候选。",
        ),
        item(
            "OVA-BODY-C4-SS3-AI03",
            "纯本体稳健备选（链数不特异）",
            "AF3 成功/失败状态差分第三二硫键 + ProteinMPNN 非共价锁环",
            "386 aa OVA P3-13R 本体点突变；无 linker、binder 或外接四聚螺旋",
            "AF3 四链同源四聚候选",
            body_ss3_aa,
            body_ss3_dna + "TAA",
            4,
            "八 seed 40/40 高置信四链连通且无 clash，ipTM/最弱链 0.77/0.76、界面 PAE 1.905 Å；29/40 完整 12 键闭环，0 非预期 S–S。C2/C3/C5/C6 ipTM 为 0.56/0.79/0.37/0.30，三聚倾向明显。",
        ),
        item(
            "OVA-BODY-C4-P208F-AI01",
            "纯本体机制备选（seed/C3 敏感）",
            "真实失败态 ProteinMPNN 负设计，以 P208F 非共价锁替代第三组 Cys",
            "386 aa OVA P3-13R 本体点突变；无 linker、binder 或外接四聚螺旋",
            "AF3 四链同源四聚候选",
            body_p208f_aa,
            body_p208f_dna + "TAA",
            4,
            "八 seed 仅 30/40 高置信且化学干净，seed 6/8 均失败；总体 ipTM/最弱链 0.795/0.80。C3 达 0.84/0.83，高于 C4，不作为首选。",
        ),
        item(
            "AI-C4-04",
            "高分但seed敏感备选",
            "ProteinMPNN C4 正选择/C3 负选择得到的四聚模块",
            "OVA P3-13R—GGGGS—NTNVILNRLTTIEAQLVRIELTLARIERLLGVK",
            "AF3 四链同源四聚候选",
            p3_aa + linker5_aa + ai_c4_04_aa,
            p3_dna + linker5_dna + ai_c4_04_dna + "TAA",
            4,
            "八 seed 共 40 个模型总体 ipTM/最弱链 0.49/0.45，仅 20/40 过阈值；虽 C2/C3 仅 0.18/0.27，但 C4 存在明显 seed 敏感性。",
        ),
        item(
            "AI-C4-06",
            "工程阳性对照（外挂螺旋）",
            "验证外接四聚螺旋可把四份 OVA 带入四聚装配",
            "OVA P3-13R—GGGGS—PTNVITARLTTIEYKLDVIELKLDRIEKLLGVK",
            "外接模块驱动的 AF3 四聚工程对照，不代表 OVA 本体界面四聚",
            p3_aa + linker5_aa + ai_c4_06_aa,
            p3_dna + linker5_dna + ai_c4_06_dna + "TAA",
            4,
            "八 seed 40 模型：ipTM/最弱链中位 0.565/0.53，35/40 过阈值、40/40 连通、0 clash；C2/C3/C5/C6 为 0.19/0.32/0.29/0.28，C4 是唯一过门槛链数。seed 8失败，仍需实验复核。",
        ),
        item(
            "CTRL-P3-pLI",
            "工程对照",
            "天然 pLI 四聚螺旋融合对照",
            "OVA P3-13R—GGGGS—GCN4-pLI",
            "四链倾向但跨 seed 不稳，仅作对照",
            p3_aa + linker5_aa + pli_aa,
            p3_dna + linker5_dna + pli_dna + "TAA",
            4,
            "首轮 C4 ipTM 0.55，但独立五-seed复核降至 0.26，且 C3 达 0.46；不作为最终候选。",
        ),
    ]

    modules = {
        "OVA WT": wt_dna,
        "OVA P3-13R": p3_dna,
        "GGGGS linker": linker5_dna,
        "GCN4-pLI": pli_dna,
        "AI-C4-04 module": ai_c4_04_dna,
        "AI-C4-06 module": ai_c4_06_dna,
        "BindCraft binder (validated mature peptide)": bindcraft_dna,
        "Stop": "TAA",
    }
    qc = [
        ("源 OVA 长度", "WT 与 P3-13R 均为 1158 bp / 386 aa，源序列无 stop。"),
        ("P3-13R 突变核验", "仅见 7 个预期氨基酸变化；从 ATG 起为 D14R/S69A/S77A/N79A/D191A/E192F/D193A。"),
        ("linker 核验", f"{len(linker5_dna)} bp / {len(linker5_aa)} aa，翻译为 {linker5_aa}。"),
        ("pLI 核验", f"{len(pli_dna)} bp / {len(pli_aa)} aa，翻译为 {pli_aa}。"),
        ("读框与终止", "所有输出 CDS 均从 ATG 开始、长度为 3 的倍数、无内部 stop，并以 TAA 终止。"),
        ("原序列保真", "WT/P3 对照及外接模块候选逐碱基保留源 OVA DNA；纯本体候选仅在各自列出的设计位点替换密码子，其余源 DNA 不变。"),
        ("GCN4-pLI 结构依据", "PDB 1GCL；平行 C4 同源四聚体。https://www.rcsb.org/structure/1GCL"),
        ("原始二聚体 AF3", "P3-13R 25 个二聚模型 ipTM 中位 0.17、界面 PAE 22.61 Å，未通过。"),
        ("OVA 实验结构对齐", "与 1.95 Å 晶体结构 1OVA 对齐：原始 P3-13R 的 50 条 AF3 链全链 Cα RMSD 中位 1.468 Å；纯本体首选的 60 条 AF3 链为 1.518 Å。"),
        ("AI-C4-04 AF3", "八 seed 40 模型总体 ipTM/最弱链 0.49/0.45，仅 20/40 过阈值；C2 0.18，C3 0.27。"),
        ("AI-C4-06 AF3", "八 seed 40 模型总体 ipTM/最弱链 0.565/0.53，35/40 过阈值、40/40 连通、0 clash；C2/C3/C5/C6 的 ipTM 为 0.19/0.32/0.29/0.28，C4 是唯一过门槛链数。"),
        ("OVA-BODY-C4-C3HOT-AI01 AF3", "纯 386 aa OVA 本体、无外接模块；八 seed 40/40 高置信四链连通、40/40 无 clash，ipTM/最弱链 0.79/0.78、界面 PAE 1.723 Å；39/40 完整干净 12 键网络。C2/C3/C5/C6 ipTM 为 0.74/0.75/0.36/0.34；C4 最高但 C2/C3 仍强，故不宣称链数特异。"),
        ("OVA-BODY-C4-SS3-AI03 AF3", "纯 386 aa OVA 本体、无外接模块；八 seed 40/40 高置信四链连通且无 clash，ipTM/最弱链 0.77/0.76、界面 PAE 1.905 Å；29/40 完整 12 键闭环。C3 ipTM 0.79，非四聚特异。"),
        ("OVA-BODY-C4-P208F-AI01 AF3", "纯 386 aa OVA 本体、无外接模块；八 seed 30/40 通过，ipTM/最弱链 0.795/0.80；C3 为 0.84/0.83，不作为首选。"),
        ("BindCraft AF3", "未融合 P3+binder 异源二聚 15/15 连通无 clash，ipTM 0.75、界面 PAE 2.57 Å；融合四聚方案目前失败。"),
        ("设计状态", "AF3 是结构计算筛选，不等于溶液中已实验证实；首选序列仍应以 SEC-MALS/AUC/native MS 复核四聚态。"),
    ]
    return constructs, modules, qc


def col_name(index: int) -> str:
    out = ""
    while index:
        index, rem = divmod(index - 1, 26)
        out = chr(65 + rem) + out
    return out


def add_cell(row: ET.Element, ref: str, value: object, style: int = 1) -> None:
    if isinstance(value, (int, float)):
        cell = ET.SubElement(row, qname(MAIN_NS, "c"), {"r": ref, "s": str(style)})
        ET.SubElement(cell, qname(MAIN_NS, "v")).text = str(value)
    else:
        cell = ET.SubElement(
            row, qname(MAIN_NS, "c"), {"r": ref, "s": str(style), "t": "inlineStr"}
        )
        inline = ET.SubElement(cell, qname(MAIN_NS, "is"))
        ET.SubElement(inline, qname(MAIN_NS, "t")).text = str(value)


def worksheet_xml(rows: list[list[object]], widths: list[float], row_heights: list[float]) -> bytes:
    root = ET.Element(qname(MAIN_NS, "worksheet"))
    ET.SubElement(root, qname(MAIN_NS, "dimension"), {"ref": f"A1:{col_name(len(widths))}{len(rows)}"})
    views = ET.SubElement(root, qname(MAIN_NS, "sheetViews"))
    view = ET.SubElement(views, qname(MAIN_NS, "sheetView"), {"workbookViewId": "0", "zoomScale": "70"})
    ET.SubElement(view, qname(MAIN_NS, "pane"), {
        "ySplit": "1", "topLeftCell": "A2", "activePane": "bottomLeft", "state": "frozen"
    })
    ET.SubElement(view, qname(MAIN_NS, "selection"), {
        "pane": "bottomLeft", "activeCell": "A2", "sqref": "A2"
    })
    ET.SubElement(root, qname(MAIN_NS, "sheetFormatPr"), {"defaultRowHeight": "15"})
    cols = ET.SubElement(root, qname(MAIN_NS, "cols"))
    for idx, width in enumerate(widths, 1):
        ET.SubElement(cols, qname(MAIN_NS, "col"), {
            "min": str(idx), "max": str(idx), "width": str(width), "customWidth": "1"
        })
    data = ET.SubElement(root, qname(MAIN_NS, "sheetData"))
    for row_idx, values in enumerate(rows, 1):
        attrs = {"r": str(row_idx)}
        if row_idx - 1 < len(row_heights):
            attrs.update({"ht": str(row_heights[row_idx - 1]), "customHeight": "1"})
        row = ET.SubElement(data, qname(MAIN_NS, "row"), attrs)
        for col_idx, value in enumerate(values, 1):
            add_cell(row, f"{col_name(col_idx)}{row_idx}", value, 2 if row_idx == 1 else 1)
    ET.SubElement(root, qname(MAIN_NS, "autoFilter"), {
        "ref": f"A1:{col_name(len(widths))}{len(rows)}"
    })
    ET.SubElement(root, qname(MAIN_NS, "pageMargins"), {
        "left": "0.25", "right": "0.25", "top": "0.5", "bottom": "0.5", "header": "0.2", "footer": "0.2"
    })
    ET.SubElement(root, qname(MAIN_NS, "pageSetup"), {"orientation": "landscape", "fitToWidth": "1"})
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def build_workbook(constructs: list[dict[str, object]], modules: dict[str, str], qc: list[tuple[str, str]]) -> None:
    candidate_headers = [
        "构建编号", "推荐级别", "设计目的", "蛋白架构", "预期装配", "AA长度/链",
        "CDS长度(bp,含stop)", "单链理论MW(kDa)", "预期装配MW(kDa)", "氨基酸序列",
        "CDS DNA 5′→3′(含TAA)", "关键风险/用途",
    ]
    candidate_rows: list[list[object]] = [candidate_headers]
    for c in constructs:
        candidate_rows.append([
            c["id"], c["priority"], c["purpose"], c["architecture"], c["assembly"], len(c["aa"]),
            len(c["cds"]), round(float(c["chain_mw"]), 2), round(float(c["assembly_mw"]), 2),
            c["aa"], c["cds"], c["risk"],
        ])
    sheet4 = worksheet_xml(
        candidate_rows,
        [20, 16, 29, 39, 42, 13, 19, 18, 19, 58, 72, 64],
        [42] + [150] * (len(candidate_rows) - 1),
    )

    module_rows: list[list[object]] = [["模块/QC项目", "AA序列或结果", "DNA序列 5′→3′", "AA长度", "DNA长度(bp)", "备注"]]
    for name, dna in modules.items():
        if name == "Stop":
            aa = "*"
        else:
            aa = translate(dna)
        note = {
            "OVA WT": "逐碱基来自源 Excel；无 stop。",
            "OVA P3-13R": "逐碱基来自源 Excel；无 stop。",
            "GGGGS linker": "与 AF3/ProteinMPNN 主候选完全一致的 5 aa C 端 linker。",
            "GCN4-pLI": "AA 采用 PDB 1GCL 序列；DNA 为通用哺乳动物偏好反向翻译。",
            "AI-C4-04 module": "ProteinMPNN C4 正选择/C3 负选择；八 seed 显示明显 seed 敏感性。",
            "AI-C4-06 module": "ProteinMPNN 当前首选模块；八 seed 中 7 个通过、35/40 模型过阈值。",
            "BindCraft binder (validated mature peptide)": "成熟肽序列；全长 P3+binder AF3 已通过，但其融合四聚体未通过。",
            "Stop": "所有输出完整 CDS 使用 TAA；如需下游标签，克隆时删除该 stop。",
        }[name]
        module_rows.append([name, aa, dna, 0 if aa == "*" else len(aa), len(dna), note])
    module_rows.append(["QC结果", "", "", "", "", "以下均由脚本自动核验"])
    for key, result in qc:
        module_rows.append([key, result, "", "", "", ""])
    sheet5 = worksheet_xml(
        module_rows,
        [25, 78, 78, 12, 15, 65],
        [42] + [115] * len(modules) + [30] + [55] * len(qc),
    )

    with ZipFile(SOURCE_XLSX) as zin:
        members = {name: zin.read(name) for name in zin.namelist()}

    workbook_root = ET.fromstring(members["xl/workbook.xml"])
    sheets_node = workbook_root.find(qname(MAIN_NS, "sheets"))
    assert sheets_node is not None
    ET.SubElement(sheets_node, qname(MAIN_NS, "sheet"), {
        "name": "四聚体候选", "sheetId": "4", qname(REL_NS, "id"): "rId7"
    })
    ET.SubElement(sheets_node, qname(MAIN_NS, "sheet"), {
        "name": "模块与QC", "sheetId": "5", qname(REL_NS, "id"): "rId8"
    })
    members["xl/workbook.xml"] = ET.tostring(workbook_root, encoding="utf-8", xml_declaration=True)

    rel_root = ET.fromstring(members["xl/_rels/workbook.xml.rels"])
    ET.SubElement(rel_root, qname(PKG_REL_NS, "Relationship"), {
        "Id": "rId7", "Type": f"{REL_NS}/worksheet", "Target": "worksheets/sheet4.xml"
    })
    ET.SubElement(rel_root, qname(PKG_REL_NS, "Relationship"), {
        "Id": "rId8", "Type": f"{REL_NS}/worksheet", "Target": "worksheets/sheet5.xml"
    })
    ET.register_namespace("", PKG_REL_NS)
    members["xl/_rels/workbook.xml.rels"] = ET.tostring(rel_root, encoding="utf-8", xml_declaration=True)

    content_root = ET.fromstring(members["[Content_Types].xml"])
    ET.SubElement(content_root, qname(CT_NS, "Override"), {
        "PartName": "/xl/worksheets/sheet4.xml",
        "ContentType": "application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml",
    })
    ET.SubElement(content_root, qname(CT_NS, "Override"), {
        "PartName": "/xl/worksheets/sheet5.xml",
        "ContentType": "application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml",
    })
    ET.register_namespace("", CT_NS)
    members["[Content_Types].xml"] = ET.tostring(content_root, encoding="utf-8", xml_declaration=True)

    app_root = ET.fromstring(members["docProps/app.xml"])
    count_node = app_root.find(f".//{qname(VT_NS, 'i4')}")
    titles_vector = app_root.find(f".//{qname(APP_NS, 'TitlesOfParts')}/{qname(VT_NS, 'vector')}")
    if count_node is not None:
        count_node.text = "5"
    if titles_vector is not None:
        titles_vector.attrib["size"] = "5"
        for title in ("四聚体候选", "模块与QC"):
            ET.SubElement(titles_vector, qname(VT_NS, "lpstr")).text = title
    ET.register_namespace("", APP_NS)
    ET.register_namespace("vt", VT_NS)
    members["docProps/app.xml"] = ET.tostring(app_root, encoding="utf-8", xml_declaration=True)

    members["xl/worksheets/sheet4.xml"] = sheet4
    members["xl/worksheets/sheet5.xml"] = sheet5
    with ZipFile(OUTPUT_XLSX, "w", ZIP_DEFLATED) as zout:
        for name, payload in members.items():
            zout.writestr(name, payload)


def build_fastas(constructs: list[dict[str, object]]) -> None:
    aa_parts = []
    cds_parts = []
    for c in constructs:
        aa_parts.append(
            f">{c['id']} | {c['architecture']} | expected: {c['assembly']}\n{wrap_sequence(str(c['aa']))}\n"
        )
        cds_parts.append(
            f">{c['id']} | {len(c['cds'])} bp | terminal_stop=TAA\n{wrap_sequence(str(c['cds']))}\n"
        )
    OUTPUT_AA.write_text("".join(aa_parts), encoding="utf-8")
    OUTPUT_CDS.write_text("".join(cds_parts), encoding="utf-8")
    top = next(c for c in constructs if c["id"] == "OVA-BODY-C4-C3HOT-AI01")
    TOP_OUTPUT_AA.write_text(
        f">{top['id']} | 386-aa pure OVA body | no linker/module\n{wrap_sequence(str(top['aa']))}\n",
        encoding="utf-8",
    )
    TOP_OUTPUT_CDS.write_text(
        f">{top['id']} | {len(top['cds'])} bp | terminal_stop=TAA\n"
        f"{wrap_sequence(str(top['cds']))}\n",
        encoding="utf-8",
    )
    ss3 = next(c for c in constructs if c["id"] == "OVA-BODY-C4-SS3-AI03")
    SS3_OUTPUT_AA.write_text(
        f">{ss3['id']} | 386-aa pure OVA body | no linker/module\n"
        f"{wrap_sequence(str(ss3['aa']))}\n",
        encoding="utf-8",
    )
    SS3_OUTPUT_CDS.write_text(
        f">{ss3['id']} | {len(ss3['cds'])} bp | terminal_stop=TAA\n"
        f"{wrap_sequence(str(ss3['cds']))}\n",
        encoding="utf-8",
    )


def main() -> None:
    constructs, modules, qc = make_constructs()
    build_workbook(constructs, modules, qc)
    build_fastas(constructs)
    for c in constructs:
        print(
            f"{c['id']}: {len(c['aa'])} aa; {len(c['cds'])} bp; "
            f"chain {c['chain_mw']:.2f} kDa; expected assembly {c['assembly_mw']:.2f} kDa"
        )
    print(f"Wrote {OUTPUT_XLSX.name}")
    print(f"Wrote {OUTPUT_AA.name}")
    print(f"Wrote {OUTPUT_CDS.name}")
    print(f"Wrote {TOP_OUTPUT_AA.name}")
    print(f"Wrote {TOP_OUTPUT_CDS.name}")
    print(f"Wrote {SS3_OUTPUT_AA.name}")
    print(f"Wrote {SS3_OUTPUT_CDS.name}")


if __name__ == "__main__":
    main()
