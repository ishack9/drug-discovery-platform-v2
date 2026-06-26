"""
Ilac Kesif Platformu v2 — admet-ai Istemcisi
admet-ai ile 50+ ADMET ozelligini lokal olarak tahmin eder.
EGFR/kanser odakli skorlama ve trafik isigi sistemi.
"""
import math
import config

# ── Ozellik Tanimi ────────────────────────────────────────────────────────────

# Kırmızı=kotu(0-3), Sari=orta(4-6), Yesil=iyi(7-10)
TRAFFIC_THRESHOLDS = {
    "hERG":                  {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "AMES":                  {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "DILI":                  {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "ClinTox":               {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "Carcinogens_Lagunin":   {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "Skin_Reaction":         {"red": 0.4,  "yellow": 0.7,  "invert": True},
    "Bioavailability_Ma":    {"red": 0.4,  "yellow": 0.6,  "invert": False},
    "HIA_Hou":               {"red": 0.4,  "yellow": 0.7,  "invert": False},
    "BBB_Martins":           {"red": 0.3,  "yellow": 0.6,  "invert": False},
    "Lipinski":              {"red": 0.4,  "yellow": 0.7,  "invert": False},
    "QED":                   {"red": 0.3,  "yellow": 0.5,  "invert": False},
    "PAMPA_NCATS":           {"red": 0.3,  "yellow": 0.6,  "invert": False},
    "Pgp_Broccatelli":       {"red": 0.4,  "yellow": 0.7,  "invert": True},
}

# EGFR/kanser icin kritik ozellikler ve agirliklari
SCORING_WEIGHTS = {
    # Toksisite (negatif — dusuk olmali)
    "hERG":                  (-2.5, True),   # Kardiyotoksisite — kritik
    "AMES":                  (-2.0, True),   # Mutajenite
    "DILI":                  (-3.0, True),   # Karaciger toksisitesi — kritik
    "ClinTox":               (-1.5, True),   # Klinik toksisite
    "Carcinogens_Lagunin":   (-1.0, True),   # Karsinojenite
    # BiyoyararlIlik (pozitif — yuksek olmali)
    "Bioavailability_Ma":    (2.0,  False),  # Oral biyoyararlilik
    "HIA_Hou":               (1.5,  False),  # Intestinal absorbsiyon
    "Lipinski":              (1.5,  False),  # Lipinski kurali
    "QED":                   (1.0,  False),  # Drug-likeness
    "PAMPA_NCATS":           (1.0,  False),  # Membran gecirgenlik
}


def _traffic_light(prop: str, value: float) -> str:
    """Ozellik degeri icin trafik isigi rengi dondurur."""
    if prop not in TRAFFIC_THRESHOLDS:
        return "gray"
    t = TRAFFIC_THRESHOLDS[prop]
    invert = t.get("invert", False)

    if invert:
        # Dusuk = iyi (toksisite icin)
        if value < t["red"]:
            return "green"
        elif value < t["yellow"]:
            return "yellow"
        else:
            return "red"
    else:
        # Yuksek = iyi (biyoyararlilik icin)
        if value > t["yellow"]:
            return "green"
        elif value > t["red"]:
            return "yellow"
        else:
            return "red"


def _compute_admet_score(props: dict) -> float:
    """
    EGFR/kanser odakli ADMET skoru hesaplar (0-10).
    Toksisite penalize eder, biyoyararlilik odul verir.
    """
    score = 7.0  # baslangic skoru

    for prop, (weight, invert) in SCORING_WEIGHTS.items():
        if prop not in props:
            continue
        val = float(props[prop])
        if math.isnan(val):
            continue

        if invert:
            # Toksisite: yuksek deger = ceza
            if val > 0.7:
                score += weight  # weight negatif
            elif val > 0.4:
                score += weight * 0.5
        else:
            # Biyoyararlilik: yuksek deger = odul
            if val > 0.7:
                score += weight * 0.5
            elif val < 0.3:
                score += weight * (-0.5)

    return round(max(0.0, min(10.0, score)), 2)


def _format_props(props: dict) -> dict:
    """Onemli ozellikleri temiz formatta dondurur."""
    key_props = {
        # Toksisite
        "hERG_kardiyotoksisite":    props.get("hERG", None),
        "AMES_mutajenite":          props.get("AMES", None),
        "DILI_karaciger":           props.get("DILI", None),
        "ClinTox_klinik":           props.get("ClinTox", None),
        "Karsino jen":              props.get("Carcinogens_Lagunin", None),
        # Biyoyararlilik
        "Oral_biyoyararlilik":      props.get("Bioavailability_Ma", None),
        "Intestinal_absorpsiyon":   props.get("HIA_Hou", None),
        "Membran_gecirgenlik":      props.get("PAMPA_NCATS", None),
        "BBB_gecis":                props.get("BBB_Martins", None),
        "Pgp_substrat":             props.get("Pgp_Broccatelli", None),
        # Drug-likeness
        "Lipinski_kurali":          props.get("Lipinski", None),
        "QED_drug_likeness":        props.get("QED", None),
        "Caco2_gecirgenlik":        props.get("Caco2_Wang", None),
        # Fizikokimyasal
        "Suda_cozunurluk":          props.get("Solubility_AqSolDB", None),
        "LogP":                     props.get("logP", None),
        "TPSA":                     props.get("tpsa", None),
        "MW":                       props.get("molecular_weight", None),
        # Metabolizma
        "CYP3A4_inhibitor":         props.get("CYP3A4_Veith", None),
        "CYP2D6_inhibitor":         props.get("CYP2D6_Veith", None),
        "Yarim_omur":               props.get("Half_Life_Obach", None),
        "LD50":                     props.get("LD50_Zhu", None),
    }
    return {k: round(float(v), 3) if v is not None and not math.isnan(float(v)) else None
            for k, v in key_props.items()}


def _traffic_summary(props: dict) -> dict:
    """Her kritik ozellik icin trafik isigi rengi dondurur."""
    raw_keys = {
        "hERG_kardiyotoksisite":  "hERG",
        "AMES_mutajenite":        "AMES",
        "DILI_karaciger":         "DILI",
        "ClinTox_klinik":         "ClinTox",
        "Oral_biyoyararlilik":    "Bioavailability_Ma",
        "Intestinal_absorpsiyon": "HIA_Hou",
        "Membran_gecirgenlik":    "PAMPA_NCATS",
        "Lipinski_kurali":        "Lipinski",
        "QED_drug_likeness":      "QED",
    }
    result = {}
    for display_key, raw_key in raw_keys.items():
        val = props.get(raw_key)
        if val is not None and not math.isnan(float(val)):
            result[display_key] = _traffic_light(raw_key, float(val))
        else:
            result[display_key] = "gray"
    return result


def evaluate_all(candidates: list) -> list:
    """
    Tum adaylari admet-ai ile degerlendirir.

    Args:
        candidates: smiles ve id iceren sozluk listesi

    Returns:
        Her aday icin admet_score, admet_props, admet_traffic eklenmis liste
    """
    try:
        from admet_ai import ADMETModel
    except ImportError:
        print("  HATA: admet-ai kurulu degil! pip install admet-ai")
        return [{**c, "admet_score": 5.0, "admet_source": "yok",
                 "admet_props": {}, "admet_traffic": {}} for c in candidates]

    print(f"  admet-ai modeli yukleniyor...")
    model = ADMETModel()

    smiles_list = [c["smiles"] for c in candidates]
    total = len(smiles_list)
    batch_size = config.ADMETLAB_BATCH_SIZE
    all_results = []

    print(f"  {total} aday tahmin ediliyor ({(total + batch_size - 1) // batch_size} batch)...")

    for i in range(0, total, batch_size):
        batch_smiles = smiles_list[i:i + batch_size]
        batch_no = i // batch_size + 1
        total_b = (total + batch_size - 1) // batch_size
        print(f"  Batch {batch_no}/{total_b} ({len(batch_smiles)} aday)...")

        try:
            df = model.predict(smiles=batch_smiles)
            for j, row in df.iterrows():
                props = row.to_dict()
                all_results.append(props)
        except Exception as e:
            print(f"  Batch {batch_no} hatasi: {e}")
            for _ in batch_smiles:
                all_results.append({})

    results = []
    for c, props in zip(candidates, all_results):
        admet_score   = _compute_admet_score(props)
        admet_props   = _format_props(props)
        admet_traffic = _traffic_summary(props)

        results.append({
            **c,
            "admet_score":   admet_score,
            "admet_source":  "admet-ai",
            "admet_props":   admet_props,
            "admet_traffic": admet_traffic,
        })

    green_count  = sum(1 for r in results if r["admet_score"] >= 7.0)
    yellow_count = sum(1 for r in results if 5.0 <= r["admet_score"] < 7.0)
    red_count    = sum(1 for r in results if r["admet_score"] < 5.0)
    print(f"  admet-ai tamamlandi: Yesil={green_count} Sari={yellow_count} Kirmizi={red_count}")

    return results