"""Additional grants database: NSF, USDA, GFI, New Harvest, DOD matches."""
import json
from datetime import datetime
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
OUT = PROJ / "p2_state_map/output"
OUT.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("ADDITIONAL GRANTS DATABASE")
print("=" * 60)

# ── Grant opportunities matched to project aims ──
grants = [
    {
        "agency": "NSF",
        "program": "EAGER",
        "title": "Early-concept Grants for Exploratory Research",
        "deadline": "Rolling",
        "amount": "$300,000 / 2 years",
        "match_score": 0.95,
        "fit_rationale": "Novel ML-driven QC methodology for emergent biotechnology; high-risk/high-reward; perfect for first 30-gene panel validation",
        "action_items": ["Contact NSF Program Officer (BIO/MCB or ENG/CBET)", "Draft 5-page EAGER proposal", "Include industry letter of interest"],
        "url": "https://www.nsf.gov/pubs/policydocs/pappg22_1/nsf22_1.pdf",
    },
    {
        "agency": "NSF",
        "program": "SBIR Phase I",
        "title": "Small Business Innovation Research",
        "deadline": "Annual (typically Oct-Dec)",
        "amount": "$275,000 / 6-12 months",
        "match_score": 0.88,
        "fit_rationale": "Commercialization-ready QC panel; requires small business entity (can partner with university via STTR)",
        "action_items": ["Identify SBIR-eligible partner or form LLC", "Draft Phase I commercialization plan"],
        "url": "https://seedfund.nsf.gov/",
    },
    {
        "agency": "USDA NIFA",
        "program": "AFRI",
        "title": "Agriculture and Food Research Initiative",
        "deadline": "Annual (typically March)",
        "amount": "$500,000-$1,000,000 / 3-4 years",
        "match_score": 0.85,
        "fit_rationale": "Alternative protein production; food safety/quality; panel directly addresses AFRI priority on sustainable agriculture",
        "action_items": ["Match to Foundational Knowledge of Agricultural Production priority area", "Partner with animal science department"],
        "url": "https://nifa.usda.gov/grants/funding-opportunities",
    },
    {
        "agency": "USDA NIFA",
        "program": "SBIR",
        "title": "USDA Small Business Innovation Research",
        "deadline": "Annual (typically Oct)",
        "amount": "$100,000 Phase I / $600,000 Phase II",
        "match_score": 0.82,
        "fit_rationale": "Food quality analytical tools; strong fit for Phase I proof-of-concept and Phase II commercialization",
        "action_items": ["Register in SAM.gov", "Partner with USDA ARS for technical merit"],
        "url": "https://nifa.usda.gov/grants/funding-opportunities/usda-nifa-sbir",
    },
    {
        "agency": "Good Food Institute",
        "program": "Research Grant Program",
        "title": "GFI Competitive Research Grant",
        "deadline": "Annual (typically June)",
        "amount": "$100,000-$250,000 / 1-2 years",
        "match_score": 0.92,
        "fit_rationale": "Directly funds cultivated meat science; GFI prioritizes analytical methods, cost reduction, and manufacturing scale-up",
        "action_items": ["Draft 3-page pre-proposal for 2027 cycle", "Emphasize cost reduction ($50/batch vs $250 RNA-seq)", "Include industry partner letter"],
        "url": "https://gfi.org/research-grants/",
    },
    {
        "agency": "Good Food Institute",
        "program": "Cultivated Meat Consortium",
        "title": "Industry-Academia Collaborative Projects",
        "deadline": "Rolling / Consortium-driven",
        "amount": "$50,000-$500,000 (cost-share)",
        "match_score": 0.80,
        "fit_rationale": "Consortium projects fund tools that benefit entire field; QC panel is broadly applicable",
        "action_items": ["Contact GFI Science & Technology team", "Propose as shared infrastructure tool"],
        "url": "https://gfi.org/consortiums/",
    },
    {
        "agency": "New Harvest",
        "program": "Research Fellowship",
        "title": "Graduate Student / Postdoc Fellowships",
        "deadline": "Annual (typically January)",
        "amount": "$30,000-$50,000 / year",
        "match_score": 0.75,
        "fit_rationale": "Funds personnel for cultivated meat research; can support grad student doing panel wet-lab validation",
        "action_items": ["Identify eligible student", "Draft research plan focused on cross-species validation"],
        "url": "https://new-harvest.org/research-fellowships/",
    },
    {
        "agency": "DOD / DARPA",
        "program": "TI",
        "title": "Tactical Innovation Office — Food Security",
        "deadline": "Rolling BAA",
        "amount": "$500,000-$2,000,000 / 2-3 years",
        "match_score": 0.70,
        "fit_rationale": "Military food security; deployed manufacturing; QC panel enables distributed production",
        "action_items": ["Frame around expeditionary/manufacturing-on-demand", "Partner with Natick Soldier Systems Center"],
        "url": "https://www.darpa.mil/work-with-us/opportunities",
    },
    {
        "agency": "DOD / AFRL",
        "program": "AFOSR",
        "title": "Air Force Office of Scientific Research — Bio-inspired Systems",
        "deadline": "Rolling / BAA",
        "amount": "$400,000-$800,000 / 3 years",
        "match_score": 0.65,
        "fit_rationale": "Bio-inspired manufacturing; tissue engineering for regenerative medicine; can pivot to defense applications",
        "action_items": ["Reframe for tissue engineering / regenerative medicine angle", "Contact AFRL program manager"],
        "url": "https://afresearchlab.com/technology/directorates/711hpw/",
    },
    {
        "agency": "ARPA-H",
        "program": "Open BAA",
        "title": "Advanced Research Projects Agency for Health",
        "deadline": "Rolling",
        "amount": "$1M-$5M / 3-5 years",
        "match_score": 0.60,
        "fit_rationale": "Panel has translational potential for muscle disease diagnostics; stretch to ARPA-H requires reframing",
        "action_items": ["Pivot to sarcopenia / muscle atrophy biomarker panel", "Partner with clinical collaborator"],
        "url": "https://arpa-h.gov/research/funding",
    },
    {
        "agency": "NSF",
        "program": "CBET",
        "title": "Biotechnology, Biochemical, and Biomanufacturing Engineering",
        "deadline": "Proposals accepted anytime (no deadline)",
        "amount": "$400,000-$600,000 / 3 years",
        "match_score": 0.90,
        "fit_rationale": "Biomanufacturing process analytical technology (PAT); NSF CBET explicitly funds bioprocess monitoring",
        "action_items": ["Draft CBET proposal focusing on process analytical technology (PAT)", "Include techno-economic analysis"],
        "url": "https://www.nsf.gov/funding/pgm_summ.jsp?pims_id=501037",
    },
    {
        "agency": "NSF",
        "program": "I-Corps",
        "title": "Innovation Corps — Commercialization Training",
        "deadline": "Cohort-based (3x/year)",
        "amount": "$50,000 / 7 weeks",
        "match_score": 0.85,
        "fit_rationale": "Team-based commercialization training; perfect next step after NSF EAGER; builds customer discovery for QC panel",
        "action_items": ["Apply with PI + entrepreneurial lead + industry mentor", "Prepare 3-minute pitch"],
        "url": "https://www.nsf.gov/news/special_reports/i-corps/",
    },
    {
        "agency": "Bill & Melinda Gates Foundation",
        "program": "Agricultural Development",
        "title": "Food Security & Nutrition",
        "deadline": "Rolling / LOI required",
        "amount": "$100,000-$1,000,000",
        "match_score": 0.72,
        "fit_rationale": "Alternative proteins for global nutrition; QC panel enables affordable, consistent production in LMICs",
        "action_items": ["Frame around affordable protein access in low-resource settings", "Partner with international agriculture institute"],
        "url": "https://www.gatesfoundation.org/about/grant-opportunities",
    },
]

# ── Summary ──
print(f"\n─── {len(grants)} Grant Opportunities ───")
print("\nBy match score (descending):")
for g in sorted(grants, key=lambda x: x["match_score"], reverse=True):
    print(f"  {g['match_score']:.2f} | {g['agency']:12} | {g['program']:15} | {g['amount']:25} | {g['deadline']}")

# ── Funding timeline ──
print("\n─── Recommended Submission Timeline ───")
timeline = [
    {"month": "Month 1", "action": "Submit NSF EAGER (rolling) — highest priority"},
    {"month": "Month 1-2", "action": "Contact GFI for pre-proposal feedback"},
    {"month": "Month 2-3", "action": "Draft USDA NIFA AFRI pre-proposal (if March deadline approaching)"},
    {"month": "Month 3-4", "action": "Submit GFI Research Grant pre-proposal"},
    {"month": "Month 4-5", "action": "Apply for NSF I-Corps cohort"},
    {"month": "Month 6-8", "action": "Draft NSF CBET full proposal (no deadline, work at own pace)"},
    {"month": "Month 8-10", "action": "Submit USDA NIFA SBIR Phase I (if Oct deadline)"},
    {"month": "Month 10-12", "action": "NSF SBIR Phase I (if Dec deadline)"},
    {"month": "Year 2", "action": "Leverage Phase I success for Phase II / follow-on grants"},
]
for t in timeline:
    print(f"  {t['month']:10} | {t['action']}")

# ── Total addressable funding ──
print("\n─── Total Addressable Funding (first 2 years) ───")
total = sum([300000, 275000, 750000, 100000, 175000, 200000, 50000, 500000, 600000, 500000, 50000])
print(f"  Conservative estimate: ${total:,}")

# ── Save ──
results = {
    "date": datetime.now().isoformat(),
    "grants": grants,
    "timeline": timeline,
    "total_addressable_2yr": total,
    "highest_priority": ["NSF EAGER", "GFI Research Grant", "NSF CBET", "USDA NIFA AFRI"],
    "personnel_grants": ["New Harvest Fellowship", "NSF GRFP (for grad student)"],
}
json.dump(results, open(OUT / "grants_database.json", "w"), indent=2)
print("\nSaved to grants_database.json")
print("DONE")
