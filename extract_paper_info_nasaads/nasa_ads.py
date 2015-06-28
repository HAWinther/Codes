import requests
import xml.etree.ElementTree as ET
import lxml.etree, os, urlparse
import time
import sys

#######################################
# Paper class that holds all information 
# about a given paper
#######################################

class Paper:
  def __init__(self):
    self.bibcode = ""
    self.title = ""
    self.authors = []
    self.abstract = ""
    self.cites = 0
    self.cites_bibcode = []

#######################################
# Generate url to ADS pages
# Input: bibcode and link_type
#        link_type = ABSTRACT, CITATIONS,
#        ...
# Output: url to desired page
#######################################

def generate_url(bibcode, link_type):
  url  = "http://adsabs.harvard.edu/cgi-bin/nph-data_query?"
  url += "bibcode=" + bibcode
  url += "&link_type=" + link_type
  url += "&db_key=AST"
  return url

#######################################
# Make a query to the NASA ADS database 
# Input: string of bibcodes in the 
#        format code1;code2;...;codeN
# Output: XML code from NASA ADS
#######################################

def nasa_ads_query(bibcodes):
  url = "http://adsabs.harvard.edu/cgi-bin/nph-abs_connect"
  print '* Performing Nasa ADS Query. Url = ', url

  # Data to be submitted via the ADS form
  form_data = {
      'bibcode':   '',
      'sort':      'CITATIONS',
      'data_type': 'XML',
      'submit':    'submit',
      }
  form_data['bibcode'] = bibcodes

  response = requests.post(url, data=form_data)
  print '* Done'
  return response.text.encode('utf-8').strip()

#######################################
# Get the HTML from a NASA ADS page
# and extract bibcodes
# Input: String with URL to library
# Output: String of bibcodes in format
#         code1;code2;...;codeN
#######################################

def extract_bibcodes_from_ads_page(url):
  print ''
  print '==================================='
  print 'Extract bibcodes from ', url
  print '==================================='
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
  print "Bibcodes = ", ";".join(bibcodes).strip()
  print '==================================='
  return ";".join(bibcodes).strip()

#######################################
# Goes through all bibcodes and extracts
# information about the paper and 
# stores it in our library
#######################################

def add_bibcodes_to_library(papers, bibcodes):
  ns = {'ref': 'http://ads.harvard.edu/schema/abs/1.1/abstracts'}
  print ''
  print '==================================='
  print 'Adding all papers to library'
  print '==================================='

  # Make query of bibcodes to NASA ADS
  xmlcode = nasa_ads_query(bibcodes)
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
    print '  * Adding paper: \"',p.title,"\" Cites: ", p.cites

  # Verbose
  print 'Total cites: ', totcites
  print 'Cites per paper: ', totcites/float(nrecords)

  # Write to file the XML
  # with open('tmp.xml', 'w') as f:
  #  f.write( xmlcode )

#######################################
# From a given bibcode extract
# the bibcodes of all the papers
# citing that paper
# Input: bibcode of paper A
# Output: bibcode(s) of all papers
#         citing paper A
#######################################

def get_bibcode_of_cites(bibcode):
  citeurl = generate_url(bibcode, "CITATIONS")
  print ''
  print '==================================='
  print "Get bibcodes of papers that cite ", bibcode
  print '==================================='
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
      print curbibcode
      bibcodes.append(curbibcode)
  return bibcodes

#######################################
# Make a list of collaborators
# Input:  A library of papers
# Output: A list of authors sorted by 
#         name
#######################################

def get_collaborators(papers):
  authors = []
  for p in papers:
    for author in p.authors:
      authors.append(author)
  return sorted(set(authors))


#######################################
# Extract publication info from NASA ADS 
# abstract page
# Input:  bibcode for a paper
# Output: publication info
#######################################

def extract_publication_info(bibcode):
  return None

#######################################
# Extract BibTeX from NASA ADS
# Input:  bibcode(s) for a paper. If more
#         then one separate them by ';'
# Output: string with BibTex
#
# NB: Assumes that '@' is only used
# to define the paper, like e.g. '@ARTICLE'
# If '@' is elsewhere in the BibTeX it 
# will fail
#######################################

def extract_bibtex(bibcodes):
  url = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode="+bibcodes+"&data_type=BIBTEX&db_key=AST"
  response = requests.get(url)
  html = response.text.encode('utf-8').strip()
  allbibtex = html.split('@')[1::]
  bibtex = "@" + "@".join(allbibtex)
  return bibtex

#######################################
# Make HTML from a list of papers 
# Input: a list of papers and the 
#        filename of the output
# Output: write filename to file and
#         returns a string with the 
#         HTML code
#######################################

def write_papers_as_html(papers,filename):
  print ''
  print '==================================='
  print 'Write HTML of library'
  print '==================================='

  totcites = 0
  totpapers = 0
  for p in papers:
    totpapers += 1
    totcites  += p.cites
  citeperpaper = totcites/float(totpapers)
  citeperpaper = int(100*citeperpaper)/100.0

  # Calculate how many unique collaborators we have
  authorlist = get_collaborators(papers)

  html  = "<html>\n<head>\n<title>My Citations</title>\n</head>\n<body>\n"
  html += "<p>\n"
  html += "  <h1>My Library</h1>\n"
  html += "  <p>\n"
  html +=     "Found a total of "+str(totpapers)+" papers with a total of "+str(totcites)+" citations. "
  html +=     "This corresponds to "+str(citeperpaper)+" citations per paper. "
  html +=     "We have a collaborated with a total of " +str(len(authorlist) - 1)+" people.\n"
  html += "  </p>\n"
  html += "</p>"
  html += "<hr>"

  # Loop over all papers
  for p in papers:
    html += "<p>\n"
    html += "  <h2>"+p.title+"</h2>\n"
    html += "  <p>\n"
    html += "    <b>Authors:</b> "
    for author in p.authors:
      # Print authors with first-name first
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
    html += "    <a href=\""+arxivhref+"\">arXiv Preprint</a>, <a href=\""+adshref+"\">NASA ADS Abstract</a>\n"
    html += "  </p>\n"
    html += "</p>\n"
    html += "<hr>\n"
  html += "</body>\n</html>"
  html = html.encode('utf-8')

  print "Write output to : ", filename
  with open(filename,"w") as f:
    f.write( html )
  return html

#####################################################################
#####################################################################
#####################################################################
#####################################################################

# Make a empty library
mylibrary = []

# Extract bibcodes from a personal NASA ADS library page
libraryurl = "http://adsabs.harvard.edu/cgi-bin/nph-abs_connect?library&libname=My+Papers&libid=5202f87fce"
bibcodes = extract_bibcodes_from_ads_page(libraryurl)

# Add all bibcodes found above to mylibrary
# and make query to NASA ADS to get relevant info
add_bibcodes_to_library(mylibrary, bibcodes)

# Get bibcodes of all papers that cite our papers
for paper in mylibrary:
  bibcode = paper.bibcode
  paper.cites_bibcode = get_bibcode_of_cites(bibcode)

# Output HTML code of all papers in your library
write_papers_as_html(mylibrary, "papers.html")

# Extract BibTeX from a list of bibcodes. Example:
bibcode = "2012ApJ...756..166W;2015MNRAS.449.3635H"
bibtex = extract_bibtex(bibcode)
print '\n BibTex = ', bibtex

