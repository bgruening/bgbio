#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import sys, os
import getopt
import numpy

from Bio.PDB.StructureBuilder import *
import xml.etree.cElementTree as etree
#import lxml.etree as etree
import warnings



class PDBMLParser:
    def __init__(self, PERMISSIVE=1, get_header=0, structure_builder=None):
        if structure_builder!=None:
            self.structure_builder=structure_builder
        else:
            self.structure_builder=StructureBuilder()
        self.header=None
        self.trailer=None
        self.line_counter=0
        self.PERMISSIVE=PERMISSIVE

        self.structure_reference = {}
        self.authors = {}
        self.current_model_id = 0
        self.current_chain_id = 0
        self.current_residue_id = 0
        self.current_atom_id = 0

        self.header_dict = {
                        'name':"", 
                        'head':'', 
                        'deposition_date' : "1900-01-01", 
                        'release_date' : "1909-01-08", 
                        'structure_method' : "unknown", 
                        'resolution' : None, 
                        'structure_reference' : "unknown", 
                        'journal_reference' : "unknown",
                        'journal': "unknown",
                        'keywords': "",
                        'author' : "", #andere Anordnung im als im normalen PDB-File, gruening, b.a. <-> b.a.gruening
                        'compound':{'1':{'misc':''}},
                        'source':{'1':{'misc':''}},
                        }


    def get_structure(self, id, file):
        """Return the structure.

        Arguments:
        o id - string, the id that will be used for the structure
        o file - name of the PDB file OR an open filehandle
        """
        self.header=None
        self.trailer=None
        # Make a StructureBuilder instance (pass id of structure as parameter)
        self.structure_builder.init_structure(id)

        self._parse(file)
        self.structure_builder.set_header(self.header)
        # Return the Structure instance
        return self.structure_builder.get_structure()

    def get_header(self):
        """Return the header."""
        return self.header

    def _parse(self, filepath):

        # get an iterable
        context = etree.iterparse(filepath, events=("start", "end"))
        structure_builder = self.structure_builder
        # turn it into an iterator
        context = iter(context)

        # get the root element
        event, root = context.next()

        header_only = False
        # TODO: its not a good solution :(
        # Save the structure state for each atom, its used to add anisotrop properties afterwards
        atom_id_structure_mapping = {}

        for event, elem in context:
            if event == "end":# and elem.tag == "record":
                if elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}atom_siteCategory":
                    elem.tail = None
                    elem.clear()

                if elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}atom_site_anisotrop":
                    #res = chain.__getitem__(("", id, ""))
                    #atom = res.__getitem__()
                    pass
                if elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}atom_site":
                    if header_only:
                        elem.clear()
                    else:
                        atom_id = int(elem.get("id"))
                        # SET the AtomID as line_counter
                        structure_builder.set_line_counter(atom_id)
                        resname = elem.find("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}auth_comp_id").text

                        if elem.find("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}group_PDB").text == "ATOM":
                            if resname=="HOH" or resname=="WAT":
                                hetero_flag="W"
                            else: 
                                hetero_flag="H"
                        else:
                            hetero_flag = " "

                        sequence_identifier = elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}auth_seq_id")
                        model = int(elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}pdbx_PDB_model_num"))
                        if self.current_model_id != model:
                            self.current_model_id = model
                            structure_builder.init_model(self.current_model_id)
                            self.current_chain_id = 0
                        seg_id = elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}label_seq_id")
                        structure_builder.init_seg(seg_id)

                        chain = elem.find("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}auth_asym_id").text
                        if self.current_chain_id != chain:
                            self.current_chain_id = chain
                            structure_builder.init_chain(self.current_chain_id)
                            self.current_residue_id = 0
                        """ 
                            Initiate a new Residue object. 

                            Arguments: 
                            o resname - string, e.g. "ASN" 
                            o field - hetero flag, "W" for waters, "H" for  
                                hetero residues, otherwise blank. 
                            o resseq - int, sequence identifier
                            o icode - string, insertion code 
                        """
                        if self.current_residue_id != sequence_identifier:
                            self.current_residue_id = sequence_identifier
                            try:
                                structure_builder.init_residue(resname, hetero_flag, sequence_identifier, " ")
                            except PDBConstructionException, message:
                                self._handle_PDB_exception(message, atom_id)
                            residue_container = (resname, hetero_flag, sequence_identifier, " ")

                        """ 
                            Initiate a new Atom object. 

                            Arguments: 
                            o name - string, atom name, e.g. CA, spaces should be stripped 
                            o coord - Numeric array (Float0, size 3), atomic coordinates 
                            o b_factor - float, B factor 
                            o occupancy - float 
                            o altloc - string, alternative location specifier 
                            o fullname - string, atom name including spaces, e.g. " CA " 
                        """
                        name = elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}label_atom_id")
                        x,y,z = float(elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}Cartn_x")), float(elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}Cartn_y")), float(elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}Cartn_z"))
                        coord = numpy.array((x, y, z), 'f')
                        b_factor = float(elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}B_iso_or_equiv", "0.0"))
                        occupancy = float(elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}occupancy", "0.0"))
                        altloc = elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}label_alt_id")
                        if not altloc:
                            altloc = ""
                        fullname = name
                        #TODO: current_atom_id wird nicht gespeichert
                        if not name:
                            print "Error: no name"
                            sys.exit()
                        try:
                            structure_builder.init_atom(name, coord, b_factor, occupancy, altloc, fullname, serial_number=atom_id)
                        except PDBConstructionException, message:
                            self._handle_PDB_exception(message, atom_id)

                        atom_container = (name, coord, b_factor, occupancy, altloc, fullname, atom_id)
                        elem.clear()
                        elem.tail = None

                        atom_id_structure_mapping[atom_id] = {"model": self.current_model_id, "seg_id": seg_id, "chain_id": self.current_chain_id, "residue": residue_container, "atom": atom_container}

                        """

                        # TODO:
                            recordtype ANISOU
                                anisou=map(float, (line[28:35], line[35:42], line[43:49], line[49:56], line[56:63], line[63:70])) 
                                # U's are scaled by 10^4 
                                anisou_array=(numpy.array(anisou, 'f')/10000.0).astype('f') 
                                structure_builder.set_anisou(anisou_array) 
                            recordtype SIGUIJ
                                # standard deviation of anisotropic B factor 
                                siguij=map(float, (line[28:35], line[35:42], line[42:49], line[49:56], line[56:63], line[63:70])) 
                                # U sigma's are scaled by 10^4 
                                siguij_array=(numpy.array(siguij, 'f')/10000.0).astype('f')    
                                structure_builder.set_siguij(siguij_array)
                            recordtype SIGATM
                                # standard deviation of atomic positions 
                                sigatm=map(float, (line[30:38], line[38:45], line[46:54], line[54:60], line[60:66]))
                                sigatm_array=numpy.array(sigatm, 'f')
                                structure_builder.set_sigatm(sigatm_array)
                        """

                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}entity_polyCategory":
                    for sub in elem.findall("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}entity_poly"):
                        self.header_dict["compound"][sub.get("entity_id", "primary")].update({"chain": sub.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}pdbx_strand_id").lower()})
                    elem.clear()
                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}entityCategory":
                    for subelem in elem.findall("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}entity"):
                        cid = subelem.get("id")
                        if self.header_dict["compound"].has_key(cid):
                            self.header_dict["compound"][cid].update({"molecule": subelem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}pdbx_description").lower()})
                        else:
                            self.header_dict["compound"][cid] = {"misc": ''}
                            self.header_dict["compound"][cid].update({"molecule": subelem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}pdbx_description").lower()})
                        temp = subelem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}pdbx_mutation")
                        if temp:
                            self.header_dict["compound"][cid].update({"mutation": temp.lower()})
                    elem.clear()
                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}entity_name_comCategory":
                    for subelem in elem:
                        cid = subelem.get("entity_id")
                        character = subelem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}name")
                        if character:
                            self.header_dict["compound"][cid].update({"synonym": character.lower()})
                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}citation_authorCategory":
                    for author in elem.findall("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}citation_author"):
                        cid = author.get("citation_id")
                        if not self.authors.has_key(cid):
                            self.authors[cid] = []
                        self.authors[cid].append(author.get("name"))
                    elem.clear()
                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}citationCategory":
                    citations = {}
                    parts = [
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}title",
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}title",
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}journal_abbrev",
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}journal_volume",
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}page_first",
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}page_last",
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}year",
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}journal_id_ISSN",
                                "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}pdbx_database_id_PubMed",
                            ]
                    for reference in elem.findall("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}citation"):
                        reference_id = reference.get("id")
                        """citation = {}
                        for child in reference.getchildren():
                            if value.character:
                                citation[key] = value.character
                        citations[reference.get("id")] = citation
                        """
                        temp_ref = []
                        for part in parts:
                            character = reference.findtext(part)
                            if character:
                                if part == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}journal_id_ISSN":
                                    temp_ref.append("issn "+character)
                                elif part == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}journal_volume":
                                    temp_ref.append("v. "+character)
                                else:
                                    temp_ref.append(character.lower())
                        self.structure_reference[reference_id] = " ".join(temp_ref)
                    #header_dict["journal_reference"] = citations
                    elem.clear()
                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}exptlCategory":
                    self.header_dict["structure_method"] = elem.find("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}exptl").get("method", "").lower()
                    elem.clear()
                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}struct_keywordsCategory":
                    for subelem in elem.getchildren():
                        self.header_dict["head"] = subelem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}pdbx_keywords", "").lower()
                        self.header_dict["keywords"] = subelem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}text", "").lower()
                    elem.clear()
                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}structCategory":
                    for subelem in elem:
                        self.header_dict["name"] = subelem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}title", "").lower()

                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}database_PDB_revCategory":
                    for rev in elem.findall("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}database_PDB_rev"):
                        if rev.get("num") == "1":
                            self.header_dict["deposition_date"] = rev.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}date_original")
                            self.header_dict["release_date"] = rev.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}date")
                    elem.clear()
                elif elem.tag == "{http://pdbml.pdb.org/schema/pdbx-v32.xsd}pdbx_entity_src_synCategory":
                    self.header_dict["source"]["1"]["organism_scientific"] = elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}organism_scientific", "").lower()
                    self.header_dict["source"]["1"]["other_details"] = elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}details", "").lower()
                    self.header_dict["source"]["1"]["organism_taxid"] = elem.findtext("{http://pdbml.pdb.org/schema/pdbx-v32.xsd}ncbi_taxonomy_id")
                    elem.clear()

                #elem.clear()
                #root.clear()
        self.header = self.header_dict
        #return structure_builder.get_structure(), header_dict

    def _handle_PDB_exception(self, message, atom_counter):
        """
        This method catches an exception that occurs in the StructureBuilder
        object (if PERMISSIVE==1), or raises it again, this time adding the 
        PDB line number to the error message.
        """
        message="%s at atom %i." % (message, atom_counter)
        if self.PERMISSIVE:
            # just print a warning - some residues/atoms will be missing
            print "PDBConstructionException: %s" % message
            print "Exception ignored.\nSome atoms or residues will be missing in the data structure."
        else:
            # exceptions are fatal - raise again with new message (including line nr)
            raise PDBConstructionException(message)



if __name__ == "__main__":
    import sys

    p=PDBMLParser(PERMISSIVE=1)

    s=p.get_structure("scr", sys.argv[1])

    for m in s:
        p=m.get_parent()
        assert(p is s)
        for c in m:
            p=c.get_parent()
            assert(p is m)
            for r in c:
                print r
                p=r.get_parent()
                assert(p is c)
                for a in r:
                    p=a.get_parent()
                    if not p is r:
                        print p, r
