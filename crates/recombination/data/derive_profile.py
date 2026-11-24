"""Development-only extraction of pinned manufacturer reference contexts; no network."""
import pathlib, hashlib, json, argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('primary_dir',type=pathlib.Path,help='Directory containing the original manufacturer sequence TXT files and cited PDF manuals/maps')
args=parser.parse_args()
root=args.primary_dir
base='https://documents.thermofisher.com/TFS-Assets/vectors/'
manual='https://assets.thermofisher.com/TFS-Assets/LSG/manuals/gatewayman.pdf'
def sha(n):return hashlib.sha256((root/n).read_bytes()).hexdigest()
def dna(n):return ''.join((root/n).read_text().split()).upper()
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
sources=[];variants=[]
for subtype,b in [(1,'ACAAGTTTGTACAAAAAAGCAGGCT'),(2,'ACCACTTTGTACAAGAAAGCTGGGT')]:
    sid=f'gateway-manual-MAN0000282-rev1.0-B{subtype}'
    sources.append(dict(id=sid,url=manual,raw_sha256=sha('gatewayman.pdf'),vector_name='Primer site string (not a vector record); manual printed page 13',vector_length=25,start_1based=1,end_1based=25,reverse_extraction=False,map_url=manual,map_sha256=sha('gatewayman.pdf')))
    variants.append(dict(id=f'attB{subtype}-primer25',family='B',subtype=subtype,sequence=b,overlap_start=8,source_id=sid))
for name,seqfile,mapfile,expected,family,coords,offset in [
    ('pDONR221','pdonr221_seq.txt','pdonr221_pdonrzeo_map.pdf',4761,'P',[(570,801,False),(2753,2984,True)],82),
    ('pDEST14','pdest14_seq.txt','pdest14_map.pdf',6422,'R',[(67,191,False),(1788,1912,True)],8),
    ('pENTR-D-TOPO','pentr_dtopo_seq.txt','pentr_dtopo_map.pdf',2580,'L',[(569,668,False),(705,804,True)],83)]:
    seq=dna(seqfile);assert len(seq)==expected and set(seq)<=set('ACGT')
    for subtype,(a,b,reverse) in enumerate(coords,1):
        sid=f'{name}-att{family}{subtype}-manufacturer-snapshot-2026-11-24'
        sources.append(dict(id=sid,url=base+seqfile,raw_sha256=sha(seqfile),vector_name=name,vector_length=len(seq),start_1based=a,end_1based=b,reverse_extraction=reverse,map_url=base+mapfile,map_sha256=sha(mapfile)))
        site=seq[a-1:b];site=rc(site) if reverse else site
        assert site[offset:offset+7]==('GTACAAA' if subtype==1 else 'GTACAAG')
        variants.append(dict(id=f'att{family}{subtype}-{name}-{len(site)}',family=family,subtype=subtype,sequence=site,overlap_start=offset,source_id=sid))
p=dict(id='manufacturer-vector-contexts-att1-att2',edition='2026-11-24.1',sources=sources,variants=variants,limitations=[
    'Explicit opt-in vector-context profile; no universal Gateway or minimal functional-site claim.',
    'Manufacturer sequence text does not provide INSDC accession/version; raw snapshot hashes and exact official map spans identify this edition.',
    'Map annotation strand is not used as matching direction: literal canonical overlap defines orientation; pENTR-D-TOPO attL1 map says complementary although reference context contains canonical GTACAAA.',
    'attR125 deletion/context variant remains distinct; LR may produce a shorter/unregistered attP hybrid; distal missing DNA is never restored.',
    'Full contexts exclude many valid vectors, distal variants and shortened sites; failure to match is unsupported applicability, not evidence of biochemical inactivity.',
    'att1.1/att2.1 and MultiSite profiles not included.'])
out=pathlib.Path(__file__).with_name('manufacturer-context-v1.json');out.write_text(json.dumps(p,indent=2)+'\n')
