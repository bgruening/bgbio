#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import sys, os
sys.path.append("..")

import unittest
from Bio.PDB.PDBParser import PDBParser
from pdbml_parser import PDBMLParser
#import xml.sax

class test_pdbml_parser(unittest.TestCase):
    def setUp(self):
        self.pdbml_structure, self.pdbml_header = pdbml_data
        self.pdb_structure, self.pdb_header = pdb_data

    def tearDown(self):
        """
            clean up used files
        """
        pass

    def test_StructureIdentity(self):
        pdb = []
        for model in self.pdb_structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        pdb.append( atom.name.strip() + str(atom.coord) + str(atom.bfactor) + str(atom.occupancy) + atom.altloc.strip() + atom.fullname.strip() )

        pdbml = []
        for model in self.pdbml_structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        pdbml.append( atom.name.strip() + str(atom.coord) + str(atom.bfactor) + str(atom.occupancy) + atom.altloc.strip() + atom.fullname.strip() )

        self.assertEqual(len(pdb), len(pdbml))
        for i in range(0,len(pdb)):
            self.assertEqual(pdb[i], pdbml[i])

        #self.assertEqual(self.pdb_structure, self.pdb_structure)#, self.pdbml_structure)

    def test_model_count(self):
        self.assertEqual(len(self.pdb_structure), len(self.pdbml_structure))

    def test_chain_count(self):
        pdb_chain_counter = 0
        for model in self.pdb_structure:
            pdb_chain_counter += len(model)
        
        pdbml_chain_counter = 0
        for model in self.pdbml_structure:
            pdbml_chain_counter += len(model)

        self.assertEqual(pdb_chain_counter, pdbml_chain_counter)


    def test_atom_count(self):
        pdb_atom_counter = 0
        for model in self.pdb_structure:
            for chain in model:
                pdb_atom_counter += len(chain)
        
        pdbml_atom_counter = 0
        for model in self.pdbml_structure:
            for chain in model:
                pdbml_atom_counter += len(chain)

        self.assertEqual(pdb_atom_counter, pdbml_atom_counter)

    def test_header(self):
        self.assertEqual(len(self.pdb_header), len(self.pdbml_header))
        for key in self.pdb_header:
            self.assertTrue(self.pdbml_header.has_key(key))
            if not key in ["author", "journal", "journal_reference"]:
                pdb_entity = self.pdb_header[key]
                if isinstance(self.pdb_header[key], basestring):
                    pdb_entity = pdb_entity.strip()
                if key in ["source", "structure_reference"]:
                    pass
                else:
                    self.assertEqual(pdb_entity, self.pdbml_header[key], "The key: '%s' is not equal! \n'%s' \n!=\n '%s'" % (key, pdb_entity, self.pdbml_header[key]))
    """
    def test_timing(self):
        import time
        print "\nPerformance ..."

        t1 = time.time()
        pdbml_file_path = "./1e74.xml"
        pdbml_file_path = "./2JOF.xml"
        a, b = pdbml_parser.run(pdbml_file_path)
        t2 = time.time()
        print 'New XML Parser took %0.3f ms' % ((t2-t1)*1000.0)

        t1 = time.time()
        pdb_file_path = "./1e74.ent"
        pdb_file_path = "./2JOF.ent"
        self.p = PDBParser()
        a,b = self.p.get_structure("a", pdb_file_path), self.p.get_header()
        t2 = time.time()
        print 'OLD PDB Parser took %0.3f ms' % ((t2-t1)*1000.0)
    """

if __name__ == "__main__":
    global pdb_data, pdbml_data


    pdb_file_path = "./1e74.ent"
    pdb_file_path = "./2JOF.ent"
    #pdb_file_path = "./data/3CMA.pdb"
    #pdb_file_path = "./data/1FVG.pdb"

    pdbml_file_path = "./1e74.xml"
    pdbml_file_path = "./2JOF.xml"
    #pdbml_file_path = "./data/3CMA.xml"
    #pdbml_file_path = "./data/1FVG.xml"

    print "PDBML-File:", pdbml_file_path
    print "PDB-File", pdb_file_path

    p = PDBMLParser()
    pdbml_data = (p.get_structure("a", pdbml_file_path), p.get_header())

    p = PDBParser()
    pdb_data = (p.get_structure("a", pdb_file_path), p.get_header())
    

    a = unittest.main()
    print "nanu\n\n---------------------------------\n\n", a
