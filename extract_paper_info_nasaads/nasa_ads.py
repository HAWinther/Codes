##############################################################################
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
##############################################################################

import requests
import xml.etree.ElementTree as ET
import lxml.etree
import time
import sys
import os.path
import urllib2

##############################################################################
# Paper class that holds all information about a given paper
##############################################################################

class Paper:
  def __init__(self):
    self.bibcode = ""
    self.title = ""
    self.authors = []
    self.abstract = ""
    self.journal = ""
    self.cites = 0
    self.cites_papers = []

##############################################################################
# Write information about a library to file (for now we do not include
# authors and abstract in file as we will mainly use it to compare cites)
# Input: list of papers and filename
# Changing this method requires a rewrite of 'check_for_new_citations()'
##############################################################################

def write_library_to_file(papers, filename, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Writing library to file: ', filename
    print '==================================================='

  data = str(len(papers)) + "\n=================================\n"
  for paper in papers:
    data += paper.bibcode + "\n" + paper.title + "\n" + str(paper.cites) + "\n"
    for p in paper.cites_papers:
      data += p.bibcode + ";"
    data += "\n=================================\n"

  with open(filename,"w") as f:
    f.write(data)

##############################################################################
# Check for new papers/cites. Read old library-file and compare with current
##############################################################################

def check_for_new_citations(papers, filename, verbose):
  oldpapers = []
  html = ""
  if(verbose):
    print ''
    print '==================================================='
    print 'Check for new papers/citations                     '
    print '==================================================='

  # Check if file exits
  if(os.path.exists(filename) == False):
    if(verbose): print 'There is no old library file (',filename,' no not exists) to check against.'
    return html

  # Read old library file and make old library
  with open(filename,"r") as f:
    data = f.readlines()
    noldpapers = int(data[0])
    for i in range(noldpapers):
      newcitebibcode = ""
      bibcode = data[2+5*i].strip()
      title = data[3+5*i].strip()
      cites = int(data[4+5*i])
      cite_bibcodes = data[5+5*i].strip().split(';')
      q = Paper()
      q.bibcode = bibcode
      q.title = title
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
  if(verbose):
    print ''
    print '==================================================='
    print 'Adding all papers to library                       '
    print '==================================================='

  ns = {'ref': 'http://ads.harvard.edu/schema/abs/1.1/abstracts'}

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

  # Verbose
  print ''
  print 'Number of papers NASA ADS =', nrecords
  print 'Total cites NASA ADS      = ', totcites
  print 'Cites per paper NASA ADS  = ', totcites/float(nrecords)

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

  # List of bibcodes is in (the first) value tag of <input type="hidden" name="bibcodes" value="...">
  for inp in root.iter("input"):
    if(inp.get("name",None)=="bibcodes" and inp.get("type",None)=="hidden"):
      bibcodes = inp.get("value",None).split(';')
      break

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
  url  = generate_url(bibcode,"ABSTRACT")
  response = requests.get(url)
  html = response.text.encode('utf-8').strip()
  root = lxml.etree.HTML(html)
  doi = arxiv = None
  for inp in root.iter("meta"):
    if(inp.get("name",None)=="dc.source"):
      journal = inp.get("content",None)
    if(inp.get("name",None)=="citation_doi"):
      doi = inp.get("content",None)
    if(inp.get("name",None)=="citation_arxiv_id"):
      arxiv = inp.get("content",None)

  pubinfo = [journal, doi, arxiv]
  return pubinfo

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
  if(totpapers>0):
    citeperpaper = totcites/float(totpapers)
    citeperpaper = int(100*citeperpaper)/100.0
  else:
    citeperpaper = 0

  # Calculate how many unique collaborators we have
  authorlist = get_collaborators(papers)
  ncollaborator = len(authorlist) - 1

  # Write HTML
  html  = "<html>\n<head>\n"
  html += "  <title>"+author+"'s Papers</title>\n"

  # Add MathJax Header
  if(useMathJax):
    html += "  <script type=\"text/x-mathjax-config\">\n MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],['\\\\(','\\\\)']]}});\n  </script>\n  <script type=\"text/javascript\" src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\"></script>\n"
  html += "</head>\n<body>\n"

  # Add introduction with some useful metrics
  html += "<h1>"+author+"'s Papers</h1>\n"
  html += "<p>\n"
  html += "  This is a list of all my papers that is listed in the NASA ADS database. We found a total of "+str(totpapers)+" papers with a total of "+str(totcites)+" citations. This corresponds to "+str(citeperpaper)+" citations per paper. The total number of collaborators on these papers is " +str(ncollaborator)+".\n"
  html += "</p>\n"
  html += "<hr>\n"

  # Add info about changes since last time we ran this code
  html += "<h2>New papers and citations:</h2>\n"
  if(newcites != ""):
    html += newcites
  else:
    html += "<p>\n No new papers/citations since last time we checked.\n</p>\n"
  html += "<hr>\n"

  # Loop over all papers in library and write info
  for p in papers:
    html += "<p>\n"

    # Add title
    html += "  <h2>"+p.title+"</h2>\n"

    # Add authors (change to firstname first)
    html += "  <p>\n"
    html += "    <b>Authors:</b> "
    for author in p.authors:
      html += author.split(',')[1] + " " + author.split(',')[0] + " ; "
    html += "\n"
    html += "  </p>\n"

    # Add abstract
    html += "  <p>\n"
    html += "    <b>Abstract:</b> "+p.abstract+"\n"
    html += "  </p>\n"

    # Add citation number
    html += "  <p>\n"
    html += "    <b>Citations (as of "+time.strftime("%d/%m/%Y")+"):</b> "+str(p.cites)+"\n"
    html += "  </p>\n"

    # Add link to journal it is published in + arXiv + Nasa ADS abstract
    pubinfo = extract_publication_info(p.bibcode)
    html   += "  <p>\n    <b>Published in: </b>"
    if(pubinfo[1] != None):
      html += "<a href=\"http://dx.doi.org/"+pubinfo[1]+"\">"+pubinfo[0]+"</a> | "
    if(pubinfo[2] != None):
      html += "<a href=\"http://arxiv.org/abs/"+pubinfo[2]+"\">eprint arXiv:"+pubinfo[2]+"</a> | "
    html   += "<a href=\""+generate_url(p.bibcode,"ABSTRACT")+"\">Nasa ADS Abstract</a>\n"
    html   += "  </p>\n"

    # Add links to paper that cite us
    if(p.cites > 0):
      html += "  <p>\n"
      html += "    <b>Links to papers that cite us: </b>"
      for pp in p.cites_papers:
        html += "<a href=\"" + generate_url(pp.bibcode,"ABSTRACT") + "\">"+pp.bibcode+"</a> | "
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

########################################
# Make Google Scholar requests. Works well,
# but as done now we only fetch the first 100 papers
########################################

def is_number(s):
  try:
    int(s)
    return True
  except ValueError:
    return False

def make_google_library(papers,username):
  url = "https://scholar.google.com/citations?user="+username
  req = urllib2.Request(url)
  response = urllib2.urlopen(req)
  html = response.read()
  root = lxml.etree.HTML(html)

  # Extract data from HTML
  npapers = 0
  hrefs   = []
  authors = [] 
  journal = []
  cites   = []
  temp    = []
  for tag in root.iter("td"):
    iclass = tag.get("class", None)
    if(iclass=="gsc_a_t"):
      for tag2 in tag.iter("a"):
        href = "https://scholar.google.com/" + tag2.get("href",None)
        hrefs.append(href)
      for tag3 in tag.iter("div"):
        temp.append(tag3.text)
    if(iclass=="gsc_a_c"):
      for tag4 in tag.iter("a"):
        if(is_number(tag4.text)): 
          cites.append(tag4.text)
  npapers = len(hrefs)
  ncites  = len(cites)

  # If 0 cites then we need to add this manually
  for i in range(npapers-ncites):
    cites.append("0")

  # Split temp array into authors and journal
  for i in range(npapers):
    authors.append(temp[2*i])
    journal.append(temp[2*i+1])

  # Make library
  totcites = 0
  for i in range(npapers):
    p = Paper()
    authorlist = authors[i].split(",")
    # If more than 8 authors google adds ... in the end
    if(len(authorlist)>= 8):
      authorlist[7] = "et al."
    p.authors = authorlist
    p.journal = journal[i]
    # We use bibcode to store the link for google
    p.bibcode = href[i]  
    p.cites = int(cites[i])
    totcites += p.cites
    papers.append(p)

  # Verbose
  print ''
  print 'Total papers in Google Scholar = ', npapers
  print 'Total cites in Google Scholar  = ', totcites
  print 'Cites per paper Google Scholar = ', totcites/float(npapers)

##############################################################################
##############################################################################
##############################################################################

# Username at google scholar
googleusername = "oDCHqpAAAAAJ"

# Create google library (Not properly tested so for now not using the library 
# for anything other than simply getting the number of papers and cites)
mygooglelibrary = []
make_google_library(mygooglelibrary,googleusername)

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
verbose = False

##############################################################################
##############################################################################
##############################################################################

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
