'''Geometry backend for ProtChemFragmentOccupancy. Reads a manifest.json and writes a CSV with the
per-molecule fragment occupancy and/or interaction scores (+ a Calpha alignment report). Needed numpy + RDKit.'''

import sys, os, csv, json
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

POLAR_CUTOFF, APOLAR_CUTOFF = 3.5, 4.0
INTERACTIONS = [('polar', {'N', 'O'}, ['N', 'O'], POLAR_CUTOFF),
                ('hydrophobic', {'C', 'S', 'CL', 'BR', 'I'}, ['C'], APOLAR_CUTOFF)]

AA3 = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
       'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'SEC', 'MSE'}

# ----------------------------------------------------------------- molecule atoms
# AutoDock atom types are not elements
AD_ELEMS = {'A': 'C', 'OA': 'O', 'NA': 'N', 'NS': 'N', 'SA': 'S', 'HD': 'H', 'HS': 'H'}

def readMol(path):
  ext = os.path.splitext(path)[1].lower()
  if ext in ('.sdf', '.mol'):
    return next((m for m in Chem.SDMolSupplier(path, removeHs=False, sanitize=False) if m is not None), None)
  if ext == '.mol2':
    # Only coordinates and elements are used: skip the substructure cleanup
    return Chem.MolFromMol2File(path, removeHs=False, sanitize=False, cleanupSubstructures=False)
  if ext == '.pdb':
    return Chem.MolFromPDBFile(path, removeHs=False, sanitize=False)

def molAtoms(path, includeH=False):
  '''(coords, elements) of the heavy (or all) atoms of the first molecule in the file, None if unreadable'''
  if path.lower().endswith('.pdbqt'):
    atoms = [(pdbElem(l, True), (float(l[30:38]), float(l[38:46]), float(l[46:54]))) for l in open(path)
             if l.startswith(('ATOM', 'HETATM'))]
  else:
    mol = readMol(path)
    if mol is None or mol.GetNumConformers() == 0:
      return None
    pos = mol.GetConformer().GetPositions()
    atoms = [(a.GetSymbol().upper(), pos[a.GetIdx()]) for a in mol.GetAtoms()]
  atoms = [(e, xyz) for e, xyz in atoms if includeH or e != 'H']
  return (np.array([xyz for _, xyz in atoms]), [e for e, _ in atoms]) if atoms else None

# ----------------------------------------------------------------- protein atoms
def pdbElem(l, pdbqt=False):
  e = l[76:79].strip().upper() or l[12:16].strip()[:1].upper()
  return AD_ELEMS.get(e, e) if pdbqt else e

def readProtein(path):
  '''Amino-acid atoms [(chain, resSeq, resName, atomName, element, xyz)] from PDB / PDBQT'''
  pdbqt = path.lower().endswith('.pdbqt')
  return [(l[21].strip() or 'A', int(l[22:26]), l[17:20].strip(), l[12:16].strip(), pdbElem(l, pdbqt),
           (float(l[30:38]), float(l[38:46]), float(l[46:54])))
          for l in open(path) if l.startswith(('ATOM', 'HETATM')) and l[17:20].strip() in AA3]

def readCA(path):
  return {(a[0], a[1]): a[5] for a in readProtein(path) if a[3] == 'CA'}

# ----------------------------------------------------------------- geometry
def kabsch(P, Q):
  '''Rotation R, translation t and RMSD mapping P onto Q (P @ R.T + t ~ Q)'''
  Pc, Qc = P - P.mean(0), Q - Q.mean(0)
  U, S, Vt = np.linalg.svd(Pc.T @ Qc)
  d = np.sign(np.linalg.det(Vt.T @ U.T))
  R = Vt.T @ np.diag([1, 1, d]) @ U.T
  rmsd = float(np.sqrt(((Qc - Pc @ R.T) ** 2).sum() / len(P)))
  return R, Q.mean(0) - P.mean(0) @ R.T, rmsd

def coverage(pose, ref, cutoff):
  '''Fraction of ref atoms with a pose atom within cutoff'''
  d2 = ((ref[:, None, :] - pose[None, :, :]) ** 2).sum(-1)
  return float((d2.min(1) <= cutoff ** 2).mean())

def interactions(lig, rec):
  '''Element-based (hydrogen-free) interaction fingerprint: {(chain, resSeq, type)} with type polar
  (N/O - N/O <= POLAR_CUTOFF, H-bond / salt bridge proxy) or hydrophobic (C/S/halogen - C <= APOLAR_CUTOFF)'''
  (xyz, elems), (recXYZ, recElems, recKeys) = lig, rec
  d = np.sqrt(((xyz[:, None, :] - recXYZ[None, :, :]) ** 2).sum(-1))
  fp = set()
  for i, e in enumerate(elems):
    for itype, ligEls, recEls, cut in INTERACTIONS:
      if e in ligEls:
        fp.update((*recKeys[j], itype) for j in np.where((d[i] <= cut) & np.isin(recElems, recEls))[0])
  return fp

def recovery(poseFP, refFP):
  '''Fraction of the reference interactions the pose reproduces'''
  return len(poseFP & refFP) / len(refFP) if refFP else 0.0

def placeFragment(fragCA, fragXYZ, recCA, recChains):
  '''One (chain, placedXYZ, rmsd, nPairs) per receptor chain, superposing the first reference chain'''
  fragChain = sorted({c for c, _ in fragCA})[0]
  out = []
  for ch in recChains:
    pairs = [(fragCA[(c, r)], recCA[(ch, r)]) for (c, r) in fragCA if c == fragChain and (ch, r) in recCA]
    if len(pairs) < 10:
      continue
    R, t, rmsd = kabsch(np.array([a for a, _ in pairs]), np.array([b for _, b in pairs]))
    out.append((ch, fragXYZ @ R.T + t, rmsd, len(pairs)))
  return out

# ----------------------------------------------------------------- main
def aggregate(agg, per, names, unionVal, threshold):
  if agg == 'union':
    return unionVal
  if agg == 'mean':
    return float(np.mean([per[n] for n in names])) if names else 0.0
  if agg == 'max':
    return max((per[n] for n in names), default=0.0)
  return float(sum(per[n] >= threshold for n in names))

def main(manifestFile):
  man = json.load(open(manifestFile))
  cutoff, agg, includeH, thr = man['cutoff'], man['aggregation'], man['includeHydrogens'], man['coverageThreshold']
  align, doInt = not man['assumeSameFrame'], man['interactions']

  rec = readProtein(man['receptor']) if align or doInt else []
  recCA = {(a[0], a[1]): a[5] for a in rec if a[3] == 'CA'}
  recChains = sorted({c for c, _ in recCA})
  if not man['placeInAllChains']:
    recChains = recChains[:1]
  if rec:
    print('Receptor: {} atoms, {} Calpha, chains {}'.format(len(rec), len(recCA), recChains))
  heavyRec = [a for a in rec if a[4] != 'H']
  recArrays = (np.array([a[5] for a in heavyRec]), np.array([a[4] for a in heavyRec]),
               [(a[0], a[1]) for a in heavyRec])

  # refs[name] = [(chain, placedCoords, interactionFP)], one per pocket
  refs, alignRows = {}, []
  for ref in man['refs']:
    name, atoms = ref['molName'], molAtoms(ref['molFile'], includeH)
    if atoms is None:
      print('WARNING: skipping reference {}: unreadable {}'.format(name, ref['molFile']))
      continue
    xyz, elems = atoms
    if not align:
      placements = [('-', xyz)]
    else:
      placed = placeFragment(readCA(ref['proteinFile']), xyz, recCA, recChains)
      if not placed:
        print('WARNING: skipping reference {}: its protein could not be aligned to the receptor'.format(name))
        continue
      for ch, _, rmsd, n in placed:
        alignRows.append((name, ch, n, round(rmsd, 3)))
        print('Reference {} -> chain {}: {} Calpha, RMSD {:.2f} A'.format(name, ch, n, rmsd))
        if n < man['minCaPairs']:
          print('WARNING: only {} matched Calpha for {} in chain {}'.format(n, name, ch))
      placements = [(ch, co) for ch, co, _, _ in placed]
    refs[name] = [(ch, co, interactions((co, elems), recArrays) if doInt else set()) for ch, co in placements]
    if doInt:
      print('Reference {}: {} interactions'.format(name, max(len(fp) for _, _, fp in refs[name])))

  refNames = list(refs)
  aggNames = [n for n in refNames if n not in set(man['exclude'])]
  # union: per pocket, the pooled heavy atoms / interactions of every aggregated reference placed there
  cloud, unionFP = {}, {}
  for n in aggNames:
    for ch, co, fp in refs[n]:
      cloud.setdefault(ch, []).append(co)
      unionFP.setdefault(ch, set()).update(fp)
  cloud = {ch: np.vstack(cos) for ch, cos in cloud.items()}

  intCols = ['fragInt'] + ['int_' + n for n in refNames] if doInt else []
  with open(man['outCsv'], 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['objId', 'molName', 'fragOcc', 'bestFrag'] + ['cov_' + n for n in refNames] + intCols)
    for pose in man['poses']:
      atoms = molAtoms(pose['poseFile'], includeH)
      if atoms is None:
        print('WARNING: unreadable pose {}'.format(pose['poseFile']))
        w.writerow([pose['objId'], pose['molName'], '', ''] + [''] * (len(refNames) + len(intCols)))
        continue

      xyz = atoms[0]
      covs = {n: max(coverage(xyz, co, cutoff) for _, co, _ in refs[n]) for n in refNames}
      unionCov = max((coverage(xyz, co, cutoff) for co in cloud.values()), default=0.0)
      best = max(aggNames, key=covs.get) if aggNames else ''
      row = [pose['objId'], pose['molName'], round(aggregate(agg, covs, aggNames, unionCov, thr), 4), best] + \
            [round(covs[n], 4) for n in refNames]

      if doInt:
        poseFP = interactions(atoms, recArrays)
        recs = {n: max(recovery(poseFP, fp) for _, _, fp in refs[n]) for n in refNames}
        unionRec = max((recovery(poseFP, fp) for fp in unionFP.values()), default=0.0)
        row += [round(aggregate(agg, recs, aggNames, unionRec, thr), 4)] + [round(recs[n], 4) for n in refNames]
      w.writerow(row)

  if align and man['alignReport']:
    with open(man['alignReport'], 'w', newline='') as f:
      w = csv.writer(f)
      w.writerow(['fragment', 'receptorChain', 'nCalpha', 'rmsd'])
      w.writerows(alignRows)


if __name__ == '__main__':
  main(sys.argv[1])
