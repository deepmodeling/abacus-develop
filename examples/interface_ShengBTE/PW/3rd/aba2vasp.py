import xml.etree.ElementTree as ET
import dpdata
import glob

def prettyXml(element, indent, newline, level = 0): 
    if element:   
        if element.text == None or element.text.isspace():     
            element.text = newline + indent * (level + 1)      
        else:    
            element.text = newline + indent * (level + 1) + element.text.strip() + newline + indent * (level + 1)    
    temp = list(element)    
    for subelement in temp:    
        if temp.index(subelement) < (len(temp) - 1):     
            subelement.tail = newline + indent * (level + 1)    
        else:  
            subelement.tail = newline + indent * level   
        prettyXml(subelement, indent, newline, level = level + 1)

def generate_xml(force, suffix):
    '''
    force : natom*3 eV/Angstrom
    '''
    root = ET.Element('modeling')
    head = ET.SubElement(root, 'calculation')
    title = ET.SubElement(head, 'varray')
    title.set('name',"forces")
    for ii in force:
        body = ET.SubElement(title, 'v')
        body.text = str(ii[0]) + ' ' + str(ii[1]) + ' ' + str(ii[2])

    tree = ET.ElementTree(root)
    prettyXml(root, '\t', '\n')

    name_xml = suffix + '/vasprun.xml'
    tree.write(name_xml, encoding = 'utf-8')

if __name__ == '__main__':
    log = glob.glob('SCF-*')
    for ii in log:
        out = dpdata.LabeledSystem(ii, fmt='abacus/scf')
        print(ii)
        generate_xml(out["forces"].reshape(-1,3), ii)