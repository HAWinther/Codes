##############################################################################
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
##############################################################################

import requests
import xml.etree.ElementTree as ET
import lxml.etree
import time
import sys
import os.path

##############################################################################
# Paper class that holds all information about a given paper
##############################################################################

class Paper:
  def __init__(self):
    self.bibcode = ""
    self.title = ""
    self.authors = []
    self.abstract = ""
    self.cites = 0
    self.cites_papers = []

def write_library_to_file(papers, filename, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Writing library to file: ', filename
    print '==================================================='

  data = str(len(papers)) + "\n"
  data   += "=================================\n"
  for paper in papers:
    data += paper.bibcode + "\n"
    data += paper.title + "\n"
    data += str(paper.cites) + "\n"
    for p in paper.cites_papers:
      data += p.bibcode + ";"
    data += "\n"
    data += "=================================\n"

  with open(filename,"w") as f:
    f.write(data)

##############################################################################
# Check for new papers/cites. Read old library-file and compare with current
##############################################################################

def check_for_new_citations(papers, filename, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Check for new papers/citations                     '
    print '==================================================='

  oldpapers = []
  html = ""

  if(os.path.exists(filename) == False):
    if(verbose): print 'There is no old library file (',filename,' no not exists)'
    return html

  # Read old library file and make old library
  with open(filename,"r") as f:
    data = f.readlines()
    noldpapers = int(data[0])
    for i in range(noldpapers):
      newcitebibcode = ""
      bibcode = data[2+5*i].strip()
      cites = int(data[4+5*i])
      cite_bibcodes = data[5+5*i].strip().split(';')
      q = Paper()
      q.bibcode = bibcode
      q.cites = cites
      for cb in cite_bibcodes:
        newpaper = Paper()
        newpaper.bibcode = cb
        q.cites_papers.append(newpaper)
      oldpapers.append(q)

  # Compare old and new library
  newpapers = len(papers) - len(oldpapers)
  newcites = 0
  if(verbose): print 'Comparing old library with current one. We have',newpapers,'new papers'

  if(newpapers>0):
    html += "<p>\n"
    html += "  We have "+str(newpapers)+" papers\n"
    html += "</p>\n"

  # Loop over all papers in the current library
  for newp in papers:
    newcitebibcode = []

    # Check if paper newp is also in old library
    match = [oldp for oldp in oldpapers if oldp.bibcode == newp.bibcode]
    if(len(match)==0):
      # Paper is not in oldpapers
      if(verbose):
        print "We found a new paper:", newp.bibcode, "that is not in the old library"
        print "This paper has",newp.cites,"citations"
      newcites += newp.cites
      for p in newp.cites_papers:
        newcitebibcode.append(p.bibcode)
    else:
      # Paper is in oldpapers. Check if citations has changed
      oldp = match[0]
      if(newp.cites > oldp.cites):
        if(verbose): print "Paper",newp.bibcode,"have",newp.cites-oldp.cites,"new citations!"
        newcites += newp.cites-oldp.cites

        # Get bibcodes of the new citations
        oldcite_bibcode = []
        for p in oldp.cites_papers:
          oldcite_bibcode.append(p.bibcode)
        newcite_bibcode = []
        for p in newp.cites_papers:
          newcite_bibcode.append(p.bibcode)
        for bibcode in newcite_bibcode:
          if bibcode not in oldcite_bibcode:
            newcitebibcode.append(bibcode)

    # Compile some info about the new cites in HTML
    if(newcitebibcode != []):
      html += "<p>\n"
      html += "  New citations for paper <b>\""+newp.title+"\"</b>: "
      for bibcode in newcitebibcode:
        adshref   = generate_url(bibcode,"ABSTRACT")
        html += "<a href=\""+adshref+"\">"+bibcode+"</a> | "
      html += "\n"
      html += "</p>\n"
  return html

##############################################################################
# Generate url to ADS pages
# Input: string with bibcode(s) [if > 1 separate by ';'] and link_type = ABSTRACT, CITATIONS, ...
# Output: url to desired page
##############################################################################

def generate_url(bibcode, link_type):
  url  = "http://adsabs.harvard.edu/cgi-bin/nph-data_query?"
  url += "bibcode="+bibcode+"&link_type="+link_type+"&db_key=AST"
  return url

def generate_bib_url(bibcode):
  url  = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?"
  url += "bibcode="+bibcode+"&data_type=BIBTEX&db_key=AST"
  return url

##############################################################################

##############################################################################
# Make a query to the NASA ADS database 
# Input: string of bibcodes in the format code1;code2;...;codeN
# Output: XML code from NASA ADS
##############################################################################

def nasa_ads_query(bibcodes, verbose):
  url = "http://adsabs.harvard.edu/cgi-bin/nph-abs_connect"
  if(verbose): print '* Performing Nasa ADS Query. Bibcodes = ', bibcodes

  # Data to be submitted via the ADS form
  form_data = {
      'bibcode':   '',
      'sort':      'CITATIONS',
      'data_type': 'XML',
      'submit':    'submit',
      }
  form_data['bibcode'] = bibcodes

  response = requests.post(url, data=form_data)
  return response.text.encode('utf-8').strip()

##############################################################################
# Get the HTML from a NASA ADS page and extract bibcodes
# Input: String with URL to library
# Output: List of bibcodes
##############################################################################

def extract_bibcodes_from_ads_page(url, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Extract bibcodes from', url
    print '==================================================='
  bibcodes = []
  response = requests.get(url)
  html = response.text.encode('utf-8').strip()
  root = lxml.etree.HTML(html)
  for inp in root.iter("input"):
    inptype = inp.get("type", None)
    inpname = inp.get("name", None)
    inpvalue = inp.get("value", None)
    if(inptype=="checkbox" and inpname=="bibcode"):
      bibcodes.append(inpvalue)
  if(verbose):
    print "Bibcodes = ", ";".join(bibcodes).strip()
    print '==================================================='
  return bibcodes

##############################################################################
# Goes through all bibcodes and extracts information about the paper and 
# stores it in our library
# Input: a (empty) list 'papers' and a list of bibcodes
##############################################################################

def add_bibcodes_to_library(papers, bibcodes, verbose):
  ns = {'ref': 'http://ads.harvard.edu/schema/abs/1.1/abstracts'}
  if(verbose):
    print ''
    print '==================================================='
    print 'Adding all papers to library                       '
    print '==================================================='

  # Make query of bibcodes to NASA ADS
  xmlcode = nasa_ads_query(bibcodes, verbose)
  root = ET.fromstring(xmlcode)
  allrecords = root.findall('ref:record',ns)
  nrecords = len(allrecords)

  # Loop over all records and extract information
  totcites = 0
  for record in allrecords:
    title    = record.find('ref:title',ns).text
    bibcode  = record.find('ref:bibcode',ns).text
    citerec  = record.find('ref:citations',ns)
    authors = []
    for author in record.findall('ref:author',ns):
      authors.append(author.text)
    cites = 0 if citerec == None else int(citerec.text)
    totcites += cites
    abstract = record.find('ref:abstract',ns).text

    # Add paper
    p = Paper()
    p.title = title
    p.authors = authors
    p.bibcode = bibcode
    p.cites = cites
    p.abstract = abstract
    papers.append(p)
    if(verbose): print '  * Adding paper: \"',p.title,"\" Cites: ", p.cites

  if(verbose):
    print 'Total cites: ', totcites
    print 'Cites per paper: ', totcites/float(nrecords)

##############################################################################
# From a given bibcode extract the bibcodes of all the papers citing that paper
# Input: bibcode of paper A
# Output: bibcode(s) of all papers citing paper A
##############################################################################

def get_bibcode_of_cites(bibcode, verbose):
  citeurl = generate_url(bibcode, "CITATIONS")
  if(verbose):
    print ''
    print '==================================================='
    print 'Get bibcodes of papers that cite', bibcode
    print '==================================================='
    print "Url = ", citeurl
  bibcodes = []
  response = requests.get(citeurl)
  html = response.text.encode('utf-8').strip()
  root = lxml.etree.HTML(html)
  for inp in root.iter("a"):
    classtype = inp.get("class", None)
    href = inp.get("href", None)
    if(classtype=="oa"):
      curbibcode = href.split('=')[1].split('&')[0]
      bibcodes.append(curbibcode)
  return bibcodes

##############################################################################
# Make a list of all collaborators we have based on all our papers
# Input:  A library of papers
# Output: A list of authors sorted by name 
##############################################################################

def get_collaborators(papers):
  authors = []
  for p in papers:
    for author in p.authors:
      authors.append(author)
  return sorted(set(authors))


##############################################################################
# Extract publication info from NASA ADS abstract page
# Input:  bibcode for a paper
# Output: string with publication info
##############################################################################

def extract_publication_info(bibcode):
  return None

##############################################################################
# Extract BibTeX from NASA ADS
# Input:  bibcode(s) for a paper. If more then one separate them by ';'
# Output: string with BibTex
#
# NB: Assumes that '@' is only used to define the paper, like e.g. '@ARTICLE'
# If '@' is elsewhere in the BibTeX it will fail!
#
# Example usage:
# bibcode = "2012ApJ...756..166W;2015MNRAS.449.3635H"
# bibtex = extract_bibtex(bibcode)
# print '\n BibTex = ', bibtex
##############################################################################

def extract_bibtex(bibcode):
  url = generate_bib_url(bibcode)
  response = requests.get(url)
  html = response.text.encode('utf-8').strip()
  allbibtex = html.split('@')[1::]
  bibtex = "@" + "@".join(allbibtex)
  return bibtex

##############################################################################
# Make HTML from a list of papers 
# Input: a list of papers and the filename of the output
# Output: write filename to file and returns a string with the HTML code
##############################################################################

def write_papers_as_html(papers, newcites, filename, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Write information about library in HTML            '
    print '==================================================='

  # Use MathJax to process math in HTML?
  useMathJax = True

  # Name of author
  author = "Hans A. Winther"

  totcites = 0
  totpapers = 0
  for p in papers:
    totpapers += 1
    totcites  += p.cites
  citeperpaper = totcites/float(totpapers)
  citeperpaper = int(100*citeperpaper)/100.0

  # Calculate how many unique collaborators we have
  authorlist = get_collaborators(papers)
  ncollaborator = len(authorlist) - 1

  html  = "<html>\n<head>\n"
  html += "<title>"+author+"'s Papers</title>\n"
  if(useMathJax):
    html += "<script type=\"text/x-mathjax-config\">\n MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],['\\\\(','\\\\)']]}});\n </script>\n <script type=\"text/javascript\" src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\"></script>\n"
  html += "</head>\n<body>\n"
  html += "<h1>"+author+"'s Papers</h1>\n"
  html += "<p>\n"
  html += "  This is a list of all my papers that is listed in the NASA ADS database. We found a total of "+str(totpapers)+" papers with a total of "+str(totcites)+" citations. This corresponds to "+str(citeperpaper)+" citations per paper. The total number of collaborators on these papers is " +str(ncollaborator)+".\n"
  html += "</p>\n"
  html += "<hr>\n"

  html += "<h2>New papers and citations:</h2>\n"
  if(newcites != ""):
    html += newcites
  else:
    html += "<p>\n No new papers/citations since last time we checked.\n</p>\n"
  html += "<hr>\n"

  # Loop over all papers
  for p in papers:
    html += "<p>\n"
    html += "  <h2>"+p.title+"</h2>\n"
    html += "  <p>\n"
    html += "    <b>Authors:</b> "
    for author in p.authors:
      # Output author-name with firstname first
      html += author.split(',')[1] + " " + author.split(',')[0] + " ; "
    html += "\n"
    html += "  </p>\n"
    html += "  <p>\n"
    html += "    <b>Abstract:</b> "+p.abstract+"\n"
    html += "  </p>\n"
    html += "  <p>\n"
    html += "    <b>Citations (as of "+time.strftime("%d/%m/%Y")+"):</b> "+str(p.cites)+"\n"
    html += "  </p>\n"
    html += "  <p>\n"
    arxivhref = generate_url(p.bibcode,"PREPRINT")
    adshref   = generate_url(p.bibcode,"ABSTRACT")
    html += "    <b>Links to paper: </b><a href=\""+arxivhref+"\">arXiv Preprint</a>, <a href=\""+adshref+"\">NASA ADS Abstract</a>\n"
    html += "  </p>\n"
    html += "  <p>\n"
    html += "    <b>Links to papers that cite us: </b>"
    for pp in p.cites_papers:
      adshref   = generate_url(pp.bibcode,"ABSTRACT")
      html += "<a href=\"" + adshref + "\">"+pp.bibcode+"</a> | "
    html += "\n"
    html += "  </p>\n"
    html += "</p>\n"
    html += "<hr>\n"
  html += "</body>\n</html>"
  html = html.encode('ascii', 'xmlcharrefreplace')

  if(verbose):
    print "Write output to : ", filename

  with open(filename,"w") as f:
    f.write( html )
  return html

##############################################################################
# Construct a library given a NASA ADS personal library. This gathers
# information about all the papers in the library and about the papers that cite it
#
# If one do not have a personal library one can replace bibcodes = ...
# by bibcodes = "bibcode1;bibcode2;...;bibcodeN" where bibcodei are valid
# NASA ADS bibcodes
##############################################################################

def construct_library(url, verbose):
  mylibrary = []
  bibcodes = ";".join( extract_bibcodes_from_ads_page(url, verbose) ).strip()
  
  # Extract all information about papers that cite us or just get bibcode?
  doquerypapersthatciteus = False

  # Add all bibcodes found above to 'papers' and make query to NASA ADS to get relevant info
  add_bibcodes_to_library(mylibrary, bibcodes, verbose)

  # Get bibcodes of all papers that cite our papers
  for paper in mylibrary:
    bibcode = paper.bibcode

    if(doquerypapersthatciteus):
      # Get all info about the papers that cite us
      bibcodes = ";".join(get_bibcode_of_cites(bibcode, verbose)).strip()
      add_bibcodes_to_library(paper.cites_papers, bibcodes, verbose)
    else:
      # Only store bibcode-data of papers that cite us in library
      for b in get_bibcode_of_cites(bibcode, verbose):
        p = Paper()
        p.bibcode = b
        paper.cites_papers.append(p)
  return mylibrary

##############################################################################
##############################################################################
##############################################################################
##############################################################################

# Filename to store the library in
libfilename = "library.txt"

# Outputfile
outputhtml = "sample_output.html"

# URL to your personal NASA ADS library
liburl = "http://adsabs.harvard.edu/cgi-bin/nph-abs_connect?library&libname=My+Papers&libid=5202f87fce"

# Verbose?
verbose = True

# Construct our library
mylibrary = construct_library(liburl, verbose)

# Check if we have new papers/citations
newcites = check_for_new_citations(mylibrary, libfilename, verbose)

# Write the library to file
write_library_to_file(mylibrary, libfilename, verbose)

# Output HTML page with info about all the papers in your library
write_papers_as_html(mylibrary, newcites, outputhtml, verbose)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
