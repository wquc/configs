# ---------------------------------------------------------------------------------
#                               USAGE OF THIS SCRIPT
# 1. Navigate to the $HOME folder, e.g. For Windows OS, it is "C:\Users\$USERNAME"
# 2. Put "pymolrc.pml" here
# 3. Put this script in the same folder, otherwise change the path in .pml file
# ---------------------------------------------------------------------------------

from pymol import cmd, util
from PIL import Image
from glob import glob
import binascii
import shutil
import os

AMINOACIDS ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
             'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', \
             'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', \
             'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}

def presetplus(obj_name="all", chain_name="A", ligand_sele="organic", label_pocket=1, label_hbdist=0):
    """
DESCRIPTION
    
    Generate customized protein-ligand interaction view.
    
SYNTAX
    
    presetplus [ obj_name [, chain_name [, ligand_sele [, label_pocket [, label_hbdist ]]]]]

EXAMPLES

    presetplus
    presetplus C
    presetplus label_hbdist = 1
    presetplus ligand_sele = sele  # After creating a selection named "sele" by clicking the ligand
    presetplus ligand_sele = sele, label_pocket = 0
    presetplus reset
    
ARGUMENTS
        
    chain_name  = string: chain identifier of target protein. If chain_name is specified as "reset", all the auxiliary objects will be deleted. {default: A} 
        
    ligand_sele = string: resn identifier of ligand inhibitor. To specify ligand name, use the HET code of ligand in PDB file. {default: determined by chemical class "organic"}
    
    label_pocket = int: label(1) or not(0) label the pocket residues. {default: 1}

    label_hbdist = int: label(1) or not(0) label the hbond distances. {default: 0}

        
NOTES
    
    1. The pocket is defined as target protein residues and water molecules (if any) within 6 Angstrom of the ligand inhibitor.

    2. The script will keep all water molecules bridging ligand and target.
    """

    # 0. delete all auxiliary objects to reset view
    cmd.delete("target")
    cmd.delete("ligand")
    cmd.delete("pocket")
    cmd.delete("hbonds")
    cmd.delete("waters")
    cmd.enable("all")
    if chain_name == "reset":
        return
        
    # 1. protein target setup
    cmd.select("target_sele", "polymer.protein and %s and chain %s"%(obj_name, chain_name))
    cmd.create("target", "target_sele")
    cmd.delete("target_sele")

    # 2. ligand inhibitor setup
    ligand_sele_str = ligand_sele if ligand_sele != "organic" else "organic"
    cmd.select("ligand_sele", ligand_sele_str)
    cmd.create("ligand", "ligand_sele")
    cmd.delete("ligand_sele")

    # 3. pocket residues setup 
    cmd.select("pocket_sele", "byres target within 6 of ligand")
    cmd.create("pocket", "pocket_sele")
    cmd.delete("pocket_sele")
    cmd.show_as("lines", "pocket")
    cmd.label("pocket and name CA", "AMINOACIDS[resn]+resi")

    # 4. hbonds setup (waters will be kept ONLY IF bridging ligand and pocket)
    hbonds_waters_ligand = cmd.find_pairs("solvent"+" & e. o", "ligand"+" & e. n+o")
    hbonds_waters_pocket = cmd.find_pairs("solvent"+" & e. o", "pocket"+" & e. n+o")
    
    idx_waters_ligand = set([str(i[0][1]) for i in hbonds_waters_ligand])
    idx_waters_pocket = set([str(i[0][1]) for i in hbonds_waters_pocket])
    idx_waters_hbonds = idx_waters_ligand.intersection(idx_waters_pocket)

    if idx_waters_hbonds:
        cmd.select("waters_sele", "byres index %s"%('+'.join(idx_waters_hbonds)))
        cmd.create("waters", "waters_sele")
        cmd.delete("waters_sele")
        cmd.show_as("nonbonded", "waters")
        cmd.distance("hbonds", "waters", "ligand", mode=2)
        cmd.distance("hbonds", "waters", "pocket", mode=2)
    cmd.distance("hbonds", "pocket", "ligand", mode=2)

    # 5. Hide everyting except manually created auxiliary objects
    cmd.disable("all")
    cmd.enable("target")
    cmd.enable("pocket")
    cmd.enable("hbonds")
    cmd.enable("waters")
    cmd.enable("ligand")
    cmd.color("gray80", "target")
    util.cbaw("pocket") # color by element
    
    # 6. label size setup of pocket and hbonds.
    if int(label_pocket):
        cmd.set("label_size", 15, "pocket")
    else:
        cmd.set("label_size", 0, "pocket")
        
    if int(label_hbdist):
        cmd.set("label_size", 15, "hbonds")
    else:
        cmd.set("label_size", 0, "hbonds")
    
    cmd.zoom("pocket")

def makemovie(movie_name="pymolmovie", duration=15):
    """
DESCRIPTION
    
    Generate movies while rotating the system.
    
SYNTAX
    
    makemovie [movie_name [, duration]]

EXAMPLES

    makemovie cdk2_view
    makemovie cdk2_view, 30
    
ARGUMENTS

    movie_name  = string: name of the movie {default: "pymolmovie"}

    duration    = float: duration of the movie in seconds {default: 15}
        
NOTES
    
    Considering performance, this command will perform 60 rotations for 6 degree/rotation.
    """

    randstr = binascii.b2a_hex(os.urandom(15)).decode("utf-8")
    tempdir = "temp-%s"%randstr
    os.mkdir(tempdir)
    
    cmd.mset("1 x60")
    util.mroll(1,60,1)
    cmd.mclear()
    cmd.mpng("%s%s%s"%(tempdir, os.sep, movie_name))
    
    frames = []
    for each_image in glob("%s%s%s*.png"%(tempdir, os.sep, movie_name)):
        each_frame = Image.open(each_image)
        frames.append(each_frame)

    duration_per_frame_in_ms = float(duration)*1000 / 60
    
    frames[0].save(movie_name+".gif", append_images=frames[1:], save_all=True, duration=duration_per_frame_in_ms, loop=0)
    shutil.rmtree(tempdir)

def interface(obj, c1='A', c2='B', cutoff=4.5, hb=1, name=""):
    """
DESCRIPTION
    
    Generate customized protein-protein interface view.
    
SYNTAX
    
    interface object, [c1 [, c2 [, cutoff[, hb[, name]]]]]

EXAMPLES

    interface 5tbm
    interface 5tbm, hb=0
    interface 5tbm, hb=0, name=before

ARGUMENTS

    obj    = string: name of complex containing more than one protein chain

    c1     = string: identifier of the 1st chain {default: A}

    c2     = string: identifier of the 2nd chain {default: B}

    cutoff = float: distance cutoff (will be ignored when hb=1) {default: 4.5}

    hb     = int: if only calculate polar contacts {default: 1}

    name   = string: identifier of contact names {default is empty}
        
NOTES
    
    cutoff = 4.5 is the default value for PyMOL preset->interface
    """

    cmd.select("tempA", "%s & polymer.protein & chain %s"%(obj, c1))
    cmd.select("tempB", "%s & polymer.protein & chain %s"%(obj, c2))
    name = "_%s"%name if len(name) else ""
    # process requests
    if not int(hb):
        contacts_name = "all_contacts%s"%name
        contacts = cmd.find_pairs("tempA", "tempB", cutoff=float(cutoff))
        index_in_contacts = "+".join(["%s+%s"%(c[0][1], c[1][1]) for c in contacts])
        cmd.create(contacts_name, "byres index %s"%index_in_contacts)
        cmd.show_as("sticks", "all_contacts%s"%name)
    else:
        contacts_res = "polar_res%s"%name
        contacts_hb  = "polar_hb%s"%name
        contacts_grp = "polar_contacts%s"%name
        contacts = cmd.find_pairs("tempA and donor", "tempB and acceptor") + \
                   cmd.find_pairs("tempB and donor", "tempA and acceptor")
        contacts = set(contacts)
        index_in_contacts = "+".join(["%s+%s"%(c[0][1], c[1][1]) for c in contacts])
        cmd.create(contacts_res, "(tempA or tempB) and byres index %s"%index_in_contacts)
        cmd.show_as("sticks", contacts_res)
        cmd.distance(contacts_hb, "tempA", "tempB", mode=2)
        cmd.group(contacts_grp, "%s %s"%(contacts_res, contacts_hb))
    # clean up    
    cmd.delete("tempA")
    cmd.delete("tempB")

def incr_gui_width():
    width = int(cmd.get("internal_gui_width"))
    if width < 1000:
        cmd.set("internal_gui_width", width+5)

def decr_gui_width():
    width = int(cmd.get("internal_gui_width"))
    if width > 100:
        cmd.set("internal_gui_width", width-5)

# COMMANDS 
cmd.extend("presetplus", presetplus)
cmd.extend("makemovie", makemovie)
cmd.extend("interface", interface)

# HOTKEYS
cmd.set_key( 'pgup' , incr_gui_width )
cmd.set_key( 'pgdn' , decr_gui_width )
