############################################################################
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
##############################################################################

import time, sys, os.path, math
import requests, urllib2 # Fetching webpages
import xml.etree.ElementTree as ET, lxml.etree # Parsing HTML/XML

##############################################################################
# Paper class that holds all information about a given paper
##############################################################################

class Paper:
  def __init__(self, authors = [], title = "", abstract = "", journal = "", 
      cites = 0, arxivcode = "", bibcode = "", inspirecode = "", googlecode = "", doi = ""):

    # Basic paper info
    self.authors  = authors
    self.title    = title
    self.abstract = abstract
    self.journal  = journal
    self.cites    = cites
  
    # List of Papers that cite this paper
    self.cites_papers = []

    # arXiv code
    self.arxivcode = arxivcode

    # Nasa ADS bibcode
    self.bibcode   = bibcode

    # The INSPIRE/HEP abstract page
    self.inspirecode = inspirecode

    # The google scholar abstract page
    self.googlecode  = googlecode

    # DOI code of paper
    self.doi = doi

    # Temp variable...
    self.code = ""

##############################################################################
# NASA/ADS Library Class
##############################################################################

class NasaAdsLibrary:
  def __init__(self, get_personal_lib=False, nasaadslibid=None, get_by_search=False, searchphrase=None, name = "My Nasa/Ads Library"):
    self.get_personal_lib = get_personal_lib
    self.nasaadslibid = nasaadslibid
    self.get_by_search = get_by_search
    self.searchphrase = searchphrase
    self.mylibrary = []
    self.newciteinfo = ""
    self.name = name

    # Creation verbose
    print ''
    print '-----> Nasa/Ads Library: name = \"',self.name,'\"'
    if(get_by_search):
      print '       Creating it by searching database for \"',self.searchphrase,'\"'
    if(get_personal_lib):
      print '       Creating it from personal library libname=\"',self.nasaadslibid[0],'\" libid=\"',self.nasaadslibid[1],'\"'

    if(get_by_search and get_personal_lib):
      print "Error: Choose to either extract a personal library *or* exctract data from search!"
      sys.exit(0)

  def construct_library(self,verbose=False):
    if(self.get_personal_lib):
      if(self.nasaadslibid == None):
        print "Error: NASA/ADS id and lib. name not provided"
        sys.exit(0)
      liburl = "http://adsabs.harvard.edu/cgi-bin/nph-abs_connect?library&libname="+self.nasaadslibid[0]+"&libid="+self.nasaadslibid[1]
    elif(self.get_by_search):
      if(self.searchphrase == None):
        print "Error: NASA/ADS searchphrase not provided"
        sys.exit(0)
      liburl = "http://adsabs.harvard.edu/cgi-bin/basic_connect?qsearch=\""+self.searchphrase+"\"&version=1"
    else:
      print "Error: Way to extract NASA/ADS library not selected!"
      sys.exit(0)
    self.mylibrary = construct_library_NasaAds(liburl, verbose)

  def compare_to_old_library(self,libfilename, verbose=False):
    self.newciteinfo = check_for_new_citations(self.mylibrary, libfilename, "NasaAds", verbose)

  def save_library(self,libfilename, verbose=False):
    write_library_to_file(self.mylibrary, libfilename, "NasaAds", verbose)

  def output_library_html(self,filename,verbose=False):
    write_library_as_html(self.mylibrary, self.name, self.newciteinfo, filename, "NasaAds", verbose)

  # Input:  bibcode(s) for a paper. If more then one bibcode separate them by ';'
  def extract_bibtex(bibcodes, verbose = False):
    url = generate_bib_url_NasaAds(bibcodes)
    html = get_html(url)
    allbibtex = html.split('@')[1::]
    bibtex = "@" + "@".join(allbibtex)
    return bibtex

##############################################################################
# Google Scholar Library Class
##############################################################################

class GoogleScholarLibrary:
  def __init__(self, googleid = None, name = "My Google Scholar Library"):
    self.googleid = googleid
    self.mylibrary = []
    self.name = name

    # Creation verbose
    print ''
    print '-----> Google Scholar Library: name = \"',self.name,'\"'
    print '       Creating it for googleid = \"',self.googleid,'\"'

  def construct_library(self, verbose = False):
    if(self.googleid == None):
      print "Error: Google Scholar Id not provided!"
      sys.exit(0)
    construct_library_GoogleScholar(self.mylibrary, self.googleid, verbose)

  def compare_to_old_library(self, libfilename, verbose = False):
    self.newciteinfo = check_for_new_citations(self.mylibrary, libfilename, "GoogleScholar", verbose)

  def save_library(self, libfilename, verbose = False):
    write_library_to_file(self.mylibrary, libfilename, "GoogleScholar", verbose)

  def output_library_html(self, filename, verbose = False):
    write_library_as_html(self.mylibrary, self.name, self.newciteinfo, filename, "GoogleScholar", verbose)

##############################################################################
# Inspire HEP Library Class
##############################################################################

class InspireHepLibrary:
  def __init__(self, inspireid = None, inspirename = None, libtype = "citeable", name = "My Inspire HEP Library"):
    self.inspireid   = inspireid
    self.inspirename = inspirename
    self.libtype     = libtype
    self.mylibrary   = []
    self.name        = name

    if(inspireid != None):
      self.inspirename = get_username_InspireHep(inspireid, False)
    elif(inspirename != None):
      self.inspireid = get_userid_InspireHep(inspirename, False)
    else:
      print "Error: We need either a inspirename or an inspire id!"
      sys.exit(0)

    # Creation verbose
    print ''
    print '-----> Inspire Hep Library: name = \"',self.name,'\"'
    print '       Creating it for id = \"',self.inspireid,'\" inspirename = \"',self.inspirename,'\"'

  def construct_library(self,verbose = False):
    if(not(self.libtype == "citeable" or self.libtype == "published")):
      print "Error: (libtype="+self.libtype+") Choose between citeable (default) or published!"
      sys.exit(0)

    construct_library_InspireHep(self.mylibrary, self.inspireid, self.inspirename, self.libtype, verbose)
    return

  def compare_to_old_library(self, libfilename, verbose = False):
    self.newciteinfo = check_for_new_citations(self.mylibrary, libfilename, "InspireHep", verbose)

  def save_library(self, libfilename, verbose=False):
    write_library_to_file(self.mylibrary, libfilename, "InspireHep", verbose)

  def output_library_html(self, filename, verbose=False):
    write_library_as_html(self.mylibrary, self.name, self.newciteinfo, filename, "InspireHep", verbose)

##############################################################################
##############################################################################
# General methods
##############################################################################
##############################################################################

##############################################################################
# Methods for fetching HTML
##############################################################################
 
def get_html(url, form_data = None):
  if(form_data == None):
    response = requests.get(url)
    html = response.text.encode('utf-8').strip()
  else:
    response = requests.post(url, data=form_data)
    html = response.text.encode('utf-8').strip()
  return html

def get_html_urllib(url):
  req = urllib2.Request(url)
  response = urllib2.urlopen(req)
  html = response.read()
  return html

##############################################################################
# Write information about a library to file (for now we do not include
# authors and abstract in file as we will mainly use it to compare cites)
# Input: list of papers and filename
# Changing this method requires a rewrite of 'check_for_new_citations()' 
##############################################################################

def papercode_from_libtype(paper, libtype):
  if(libtype=="NasaAds"):
    code = paper.bibcode
  elif(libtype=="InspireHep"):
    code = paper.inspirecode
  elif(libtype=="GoogleScholar"):
    code = paper.googlecode
  else:
    print 'Error: unknown libtype=',libtype
    sys.exit(0)
  return code

def write_library_to_file(papers, filename, libtype, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Writing', libtype, 'library to file:', filename
    print '==================================================='

  data = str(len(papers)) + "\n=================================\n"
  for paper in papers:
    code = papercode_from_libtype(paper, libtype)
    data += code + "\n" + paper.title + "\n" + str(paper.cites) + "\n"
    for p in paper.cites_papers:
      code = papercode_from_libtype(p, libtype)
      data += code + ";"
    data += "\n=================================\n"

  with open(filename,"w") as f:
    f.write(data)

##############################################################################
# Make HTML from a list of papers
# Input: a list of papers and the filename of the output
# Output: write filename to file and returns a string with the HTML code
##############################################################################

def write_library_as_html(papers, libname, newcites, filename, libtype, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Write information about library in HTML            '
    print '==================================================='

  # Use MathJax to process math in HTML?
  useMathJax = True

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
  authorlist = get_all_collaborators(papers)
  ncollaborator = len(authorlist) - 1

  # Write HTML
  html  = "<html>\n<head>\n"
  html += "  <title>"+libname+"</title>\n"

  # Add MathJax Header
  if(useMathJax):
    html += "  <script type=\"text/x-mathjax-config\">\n MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],['\\\\(','\\\\)']]}});\n  </script>\n  <script type=\"text/javascript\" src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\"></script>\n"
  html += "</head>\n<body>\n"

  # Add introduction with some useful metrics
  html += "<h1>"+libname+"</h1>\n"
  html += "<p>\n"
  html += "  We found a total of "+str(totpapers)+" papers with a total of "+str(totcites)+" citations. This corresponds to "+str(citeperpaper)+" citations per paper. The total number of collaborators on these papers is " +str(ncollaborator)+".\n"
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
      tmp = ""
      try:
        tmp = author.split(',')[1] + " " + author.split(',')[0] + " ; "
      except:
        tmp = author + " ; "
        pass
      html += tmp
    html += "\n"
    html += "  </p>\n"

    # Add abstract
    html += "  <p>\n"
    tmp = ""
    try:
      tmp = "    <b>Abstract:</b> "+p.abstract+"\n"
    except:
      tmp = "    <b>Abstract:</b> No abstract availiable\n"
      pass
    html += tmp
    html += "  </p>\n"

    # Add citation number
    html += "  <p>\n"
    html += "    <b>Citations (as of "+time.strftime("%d/%m/%Y")+"):</b> "+str(p.cites)+"\n"
    html += "  </p>\n"

    # Add link to journal it is published in + arXiv + [Nasa ADS abstract / Inspire Hep / Google Scholar]
    html   += "  <p>\n    <b>Published in: </b>"
    if(p.doi != ""):
      html += "<a href=\"http://dx.doi.org/"+p.doi+"\">"+(p.doi if p.journal == "" else p.journal)+"</a> | "
    if(p.arxivcode != ""):
      html += "<a href=\"http://arxiv.org/abs/"+p.arxivcode+"\">eprint arXiv:"+p.arxivcode+"</a> | "
    if(p.bibcode != ""):
      html += "<a href=\""+generate_url_NasaAds(p.bibcode,"ABSTRACT")+"\">Nasa ADS Abstract</a>\n"
    if(p.inspirecode != ""):
      html += "<a href=\""+generate_url_InspireHep(p.inspirecode, 'abstract')+"\">Inspire HEP</a>\n"
    html   += "  </p>\n"

    # Add links to paper that cite us
    if(p.cites > 0):
      html += "  <p>\n"
      html += "    <b>Links to papers that cite us: </b>"
      for i,pp in enumerate(p.cites_papers):
        if(libtype=="NasaAds"):
          html += "<a href=\"" + generate_url_NasaAds(pp.bibcode,"ABSTRACT") + "\">["+str(i+1)+"]</a>"
        elif(libtype=="InspireHep"):
          html += "<a href=\"" + generate_url_InspireHep(pp.inspirecode, 'abstract') + "\">["+str(i+1)+"]</a>"
        elif(libtype=="GoogleScholar"):
          html += "<a href=\"https://scholar.google.no/citations?view_op=view_citation&citation_for_view=" + pp.googlecode + "\">["+str(i+1)+"]</a>"
        html += " | " if i<p.cites-1 else ""
      html += "\n"
      html += "  </p>\n"

    html += "</p>\n"
    html += "<hr>\n"
  html += "</body>\n</html>"

  if(libtype=="NasaAds"):
    html = html.encode('ascii', 'xmlcharrefreplace')

  if(verbose):
    print "Write output to : ", filename

  with open(filename,"w") as f:
    f.write( html )

  return html

##############################################################################
# Make a list of all collaborators we have based on all our papers
# Input:  A library of papers
# Output: A list of authors sorted by name
##############################################################################

def get_all_collaborators(papers):
  authors = []
  for p in papers:
    for author in p.authors:
      authors.append(author)
  return sorted(set(authors))

##############################################################################
##############################################################################
# Inspire HEP specific methods
##############################################################################
##############################################################################

##############################################################################
# Find the Inspire ID of a given username
##############################################################################

def get_userid_InspireHep(inspirename, verbose):
  url = "http://inspirehep.net/author/profile/"+inspirename
  html = get_html(url)
  root = lxml.etree.HTML(html)
  for tag in root.findall(".//div[@class='row'][@id='profile-page'][@data-box-pid]"):
    inspireid = tag.get("data-box-pid",None)
  if(verbose): print 'The userid for',inspirename,'is',inspireid
  return inspireid

def get_username_InspireHep(inspireid, verbose):
  url = "http://inspirehep.net/author/profile/"+inspireid
  html = get_html(url)
  root = lxml.etree.HTML(html)
  for tag in root.findall(".//li[@class='active']/a"):
    inspirename = tag.get("href",None).split("/profile/")[1]
  if(verbose): print 'The username for',inspireid,'is',inspirename
  return inspirename

##############################################################################
# Fetch data from Inspire Hep given a userid
##############################################################################

def construct_library_InspireHep(papers, inspireid, inspirename, libtype, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Construct Inspire Hep library', inspirename, inspireid, libtype
    print '==================================================='

  # Extract all information about papers that cite us or just get inspirecode?
  doquerypapersthatciteus = False

  # Fetch citation summary
  citeurl = "http://inspirehep.net/author/profile/citations-summary"
  form_data = {'jsondata': '{\"personId\":'+str(inspireid)+'}',}
  html = get_html(citeurl, form_data)

  root = lxml.etree.HTML(html)
  tags = root.findall(".//tr/td")
  cites     = int(tags[3].text)
  cites_pub = int(tags[4].text)

  # Fetch papers (Inspire only display 250 at the time)
  for i in range(1000):
    url = "http://inspirehep.net/search?p=author:"+inspirename+" AND collection:"+libtype+"&rg=250&jrec="+str(i*250+1)
    html = get_html(url)
    root = lxml.etree.HTML(html)

    if(verbose):
      print 'Fetching :', url

    # Extract paper IDs
    for tag in root.iter("abbr"):
      tagclass = tag.get("class", None)
      if(tagclass=="unapi-id"):
        p = Paper()
        p.inspirecode = tag.get("title",None)
        papers.append(p)
      
    # Fetch arXiv number if it exits
    for i,record in enumerate(root.findall(".//div[@class='record_body']")):
      for tag in record.findall(".//small/b/a"):
        if(tag.text != None):
          if(tag.text.find("arXiv:") >=0): 
            papers[i].arxivcode = tag.text.split("arXiv:")[1]

    if(i >= int( math.floor( (cites-0.5)/250.0 ) ) ): break

  # From the paper ID extract XML page with authors, title and abstract
  for paper in papers: 
    extract_main_paper_data_InspireHep(paper, paper.inspirecode, verbose)

  # Fetch citations and add to list
  for paper in papers: 
    cites_inspirecode = get_inspirecodes_of_cites(paper.inspirecode, verbose)
    paper.cites = len(cites_inspirecode)
    for c in cites_inspirecode:
      p = Paper()
      p.inspirecode = c
      paper.cites_papers.append(p)

      # Fetch more info about the paper that cite us
      if(doquerypapersthatciteus):
        extract_main_paper_data_InspireHep(p, c, False)

  # Verbose
  print ''
  print '       Total papers in Inspire HEP    = ', len(papers)
  print '       Total cites  in Inspire HEP    = ', cites,'(',cites_pub,')'
  print '       Cites per paper Inspire HEP    = ', float("{0:.2f}".format( cites/float(len(papers)) ))

##############################################################################
# Extracts title, abstract, authors and doi (if it exists) from a inspirecode
# and stores it in paper
##############################################################################

def extract_main_paper_data_InspireHep(paper, inspirecode, verbose):
    url = generate_url_InspireHep(paper.inspirecode,"xml")
    xml = get_html(url)
    root = lxml.etree.HTML(xml) 
 
    # xxx should not be needed but is!
    paper.authors = []
    paper.title = ""
    paper.abstract = ""

    if(verbose):
      print ''
      print 'Extracting data for inspirecode = ', inspirecode
      print 'Url = ', url

    for author   in root.iter("author"):
      if author.text   != None: paper.authors.append(author.text)
    for abstract in root.iter("abstract"):
      if abstract.text != None: paper.abstract = abstract.text.encode('utf-8')
    for title    in root.iter("title"):
      if title.text    != None: paper.title = title.text.encode('utf-8')
    for doi      in root.iter("electronic-resource-num"):
      if doi.text      != None: paper.doi = doi.text

    if(verbose):
      print paper.title,'|',paper.authors,'|',paper.doi

##############################################################################
# Generate different url's related to an Inspire HEP entry
##############################################################################

def generate_url_InspireHep(inspirecode, urltype = ''):
  dictionary = {'bibtex': 'export/hx', 'latex': 'export/hlxu', 'xml': 'export/xe', 'harvmac': 'hlxh', 'abstract': '', '': '', 'citations': 'citations', 'marcxml': 'export/xm', 'marc': 'export/hm', 'nml': 'xn','dc': 'xd'}

  code = urltype
  for key in dictionary:
    code = code.replace(key,dictionary[key])
  
  url = "http://inspirehep.net/record/"+inspirecode+"/"+code
  return url

##############################################################################
# Collect information about all papers that cite a given paper
##############################################################################

def get_inspirecodes_of_cites(inspirecode, verbose):
  cites_inspirecode = []
  cites = 0

  if(verbose):
    print ''
    print '==================================================='
    print 'Fetching inspirecodes of papers that cite',inspirecode
    print '==================================================='

  # Loop over different pages containing data (since 250 max per page)
  for i in range(1000):
    url = "http://inspirehep.net/search?ln=en&ln=en&p=refersto:recid:"+inspirecode+"&of=hb&action_search=Search&rg=250&jrec="+str(i*250+1)

    html = get_html(url)
    root = lxml.etree.HTML(html)

    # Get number of cites
    if i==0:
      for c in root.findall(".//td[@class='searchresultsboxheader'][@align='center']/strong"):
        cites = int(c.text)

    # Get all inspirecodes from this page
    for a in root.findall(".//div[@class='record_body']/a[@class='titlelink']"):
      cites_inspirecode.append(a.get("href",None).split("record/")[1])
  
    if(i >= int( math.floor( (cites-0.5)/250.0 ) ) ): break

  if(verbose): 
    print 'We found',len(cites_inspirecode),'inspire codes (total cites = ',cites,')'

  return cites_inspirecode

##############################################################################
##############################################################################
# Nasa ADS specific methods
##############################################################################
##############################################################################

##############################################################################
# Check for new papers/cites. Read old library-file and compare with current
##############################################################################

def check_for_new_citations(papers, filename, libtype, verbose):

  # Copy to paper-code to the 'code' variable
  for p in papers:
    p.code = papercode_from_libtype(p, libtype)
    for c in p.cites_papers:
      c.code = papercode_from_libtype(c, libtype)
  
  oldpapers = []
  html = ""
  if(verbose):
    print ''
    print '==================================================='
    print 'Check for new papers/citations [',libtype,']       '
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
      newcitecode = ""
      code = data[2+5*i].strip()
      title = data[3+5*i].strip()
      cites = int(data[4+5*i])
      cite_codes = data[5+5*i].strip().split(';')

      q = Paper(title=title,cites=cites)
      q.code = code

      for cb in cite_codes:
        newpaper = Paper()
        newpaper.code = cb
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
    newcitecode = []

    # Check if paper newp is also in old library
    match = [oldp for oldp in oldpapers if oldp.code == newp.code]

    if(len(match)==0):
      # Paper is not in oldpapers
      if(verbose):
        print "We found a new paper:", newp.code, "that is not in the old library"
        print "This paper has",newp.cites,"citations"
      newcites += newp.cites
      for p in newp.cites_papers:
        newcitecode.append(p.code)
    else:
      # Paper is in oldpapers. Check if citations has changed
      oldp = match[0]
      if(newp.cites > oldp.cites):
        if(verbose): print "Paper",newp.code,"have",newp.cites-oldp.cites,"new citations!"
        newcites += newp.cites-oldp.cites

        # Get bibcodes of the new citations
        oldcite_code = []
        for p in oldp.cites_papers:
          oldcite_code.append(p.code)
        newcite_code = []
        for p in newp.cites_papers:
          newcite_code.append(p.code)
        for code in newcite_code:
          if code not in oldcite_code:
            newcitecode.append(code)

    # Compile some info about the new cites in HTML
    if(newcitecode != []):
      html += "<p>\n"
      html += "  New citations for paper <b>\""+newp.title+"\"</b>: "
      for code in newcitecode:
        if(libtype=="NasaAds"):
          abshref = generate_url_NasaAds(code,"ABSTRACT")
        elif(libtype=="InspireHep"):
          abshref = generate_url_InspireHep(code, 'abstract')
        elif(libtype=="GoogleScholar"):
          abshref = "https://scholar.google.no/citations?view_op=view_citation&citation_for_view=" + code
        html += "<a href=\""+abshref+"\">"+code+"</a> | "
      html += "\n"
      html += "</p>\n"
  return html

##############################################################################
# Generate url to Nasa/ADS pages
# Input: string with bibcode(s) [if > 1 separate by ';'] and link_type = ABSTRACT, CITATIONS, ...
# Output: url to desired page
##############################################################################

def generate_url_NasaAds(bibcode, link_type):
  url  = "http://adsabs.harvard.edu/cgi-bin/nph-data_query?"
  url += "bibcode="+bibcode+"&link_type="+link_type+"&db_key=AST"
  return url

def generate_bib_url_NasaAds(bibcode):
  url  = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query?"
  url += "bibcode="+bibcode+"&data_type=BIBTEX&db_key=AST"
  return url

##############################################################################
# Make a query to the NASA ADS database
# Input: string of bibcodes in the format code1;code2;...;codeN
# Output: XML code from NASA ADS
##############################################################################

def query_NasaAds(bibcodes, verbose):
  url = "http://adsabs.harvard.edu/cgi-bin/nph-abs_connect"
  if(verbose): print '* Performing Nasa ADS Query. Bibcodes = ', bibcodes
  form_data = {
      'bibcode':   '',
      'sort':      'CITATIONS',
      'data_type': 'XML',
      'submit':    'submit',
      }
  form_data['bibcode'] = bibcodes
  html = get_html(url, form_data)
  return html

##############################################################################
# Get the HTML from a NASA ADS page and extract bibcodes
# Input: String with URL to library
# Output: List of bibcodes
##############################################################################

def extract_bibcodes_from_page_NasaAds(url, verbose):
  bibcodes = []
  html = get_html(url)
  root = lxml.etree.HTML(html)
  for i in root.findall(".//input[@type='checkbox'][@name='bibcode']"):
    bibcodes.append(i.get("value",None))

  if(verbose):
    print ''
    print '==================================================='
    print 'Extract bibcodes from', url
    print '==================================================='
    print "Bibcodes = ", ";".join(bibcodes).strip()

  return bibcodes

##############################################################################
# Goes through all bibcodes and extracts information about the paper and
# stores it in our library
# Input: a (empty) list 'papers' and a list of bibcodes
##############################################################################

def add_bibcodes_to_library_NasaAds(papers, bibcodes, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Adding all papers to library                       '
    print '==================================================='

  ns = {'ref': 'http://ads.harvard.edu/schema/abs/1.1/abstracts'}

  # Make query of bibcodes to NASA ADS
  xmlcode = query_NasaAds(bibcodes, verbose)
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
    p = Paper(title=title,authors=authors,bibcode=bibcode,cites=cites,abstract=abstract)
    extract_publication_info_NasaAds(p, verbose)
    papers.append(p)
    if(verbose): print '  * Adding paper: \"',p.title,"\" Cites: ", p.cites

  # Verbose
  print ''
  print '       Total papers in Nasa/ADS       = ', nrecords
  print '       Total cites  in Nasa/ADS       = ', totcites
  print '       Cites per paper Nasa/ADS       = ', float("{0:.2f}".format( totcites/float(nrecords) ))

##############################################################################
# From a given bibcode extract the bibcodes of all the papers citing that paper
# Input: bibcode of paper A
# Output: bibcode(s) of all papers citing paper A
##############################################################################

def get_bibcodes_of_cites_NasaAds(bibcode, verbose):
  citeurl = generate_url_NasaAds(bibcode, "CITATIONS")
  if(verbose):
    print ''
    print '==================================================='
    print 'Get bibcodes of papers that cite', bibcode
    print '==================================================='
    print "Url = ", citeurl

  html = get_html(citeurl)
  root = lxml.etree.HTML(html)

  bibcodes = []
  for inp in root.iter("input"):
    if(inp.get("name",None)=="bibcodes" and inp.get("type",None)=="hidden"):
      bibcodes = inp.get("value",None).split(';')
      break

  return bibcodes

##############################################################################
# Extract publication info from NASA ADS abstract page
# Input:  a paper (bibcode need to be defined)
##############################################################################

def extract_publication_info_NasaAds(paper, verbose):
  if(paper.bibcode == ""):
    print "Error in extract_publication_info_NasaAds. Bibcode not defined!"
    sys.exit(0)

  if(verbose):
    print 'extract_publication_info_NasaAds for bibcode = ', paper.bibcode

  url  = generate_url_NasaAds(paper.bibcode,"ABSTRACT")
  html = get_html(url)
  root = lxml.etree.HTML(html)
  doi = arxiv = None

  for inp in root.iter("meta"):
    if(inp.get("name",None)=="dc.source"):
      journal = inp.get("content",None)
      paper.journal = journal if journal != None else ""
    if(inp.get("name",None)=="citation_doi"):
      doi = inp.get("content",None)
      paper.doi = doi if doi != None else ""
    if(inp.get("name",None)=="citation_arxiv_id"):
      arxivcode = inp.get("content",None)
      paper.arxivcode = arxivcode if arxivcode != None else ""

##############################################################################
# Construct a library given a NASA ADS personal library. This gathers
# information about all the papers in the library and about the papers that cite it
#
# If one do not have a personal library one can replace bibcodes = ...
# by bibcodes = "bibcode1;bibcode2;...;bibcodeN" where bibcodei are valid
# NASA ADS bibcodes
##############################################################################

def construct_library_NasaAds(url, verbose):
  mylibrary = []
  bibcodes = ";".join( extract_bibcodes_from_page_NasaAds(url, verbose) ).strip()

  if(verbose):
    print ''
    print '==================================================='
    print 'Construct library Nasa ADS from a personal library'
    print '==================================================='
    print 'Fetch library url:', url

  # Extract all information about papers that cite us or just get bibcode?
  doquerypapersthatciteus = False

  # Add all bibcodes found above to 'papers' and make query to NASA ADS to get relevant info
  add_bibcodes_to_library_NasaAds(mylibrary, bibcodes, verbose)

  # Get bibcodes of all papers that cite our papers
  for paper in mylibrary:
    bibcode = paper.bibcode

    if(doquerypapersthatciteus):
      bibcodes = ";".join(get_bibcodes_of_cites_NasaAds(bibcode, verbose)).strip()
      add_bibcodes_to_library_NasaAds(paper.cites_papers, bibcodes, verbose)
    else:
      # Only store bibcode-data of papers that cite us in library
      for b in get_bibcodes_of_cites_NasaAds(bibcode, verbose):
        p = Paper(bibcode=b)
        paper.cites_papers.append(p)
  return mylibrary

##############################################################################
##############################################################################
# Google Scholar specific methods
##############################################################################
##############################################################################

##############################################################################
# Make Google Scholar requests. Works well, but as done now we only fetch 
# the first 100 papers on the page! Addition requests (form submits) are needed
# to extract the rest...
##############################################################################

def construct_library_GoogleScholar(papers, googleid, verbose):
  url = "https://scholar.google.com/citations?user="+googleid
  html = get_html_urllib(url)
  
  hrefs = []; authors = []; journal = []; cites = []; temp = []; titles = []
  
  if(verbose or True):
    print ''
    print '==================================================='
    print 'Construct Google Scholar library', googleid
    print '==================================================='
    print 'Fetch url:', url

  # Extract data from HTML
  root = lxml.etree.HTML(html)
  for tag in root.findall(".//td[@class='gsc_a_t']/a"):
    hrefs.append("https://scholar.google.com/"+tag.get("href",None))
    titles.append(tag.text)
  for tag in root.findall(".//td[@class='gsc_a_t']/div"):
    temp.append(tag.text)
  for tag in root.findall(".//td[@class='gsc_a_c']/a"):
    cites.append(tag.text if tag.text.strip() != "" else "0")
  npapers = len(hrefs)
  for i in range(npapers):
    authors.append(temp[2*i])
    journal.append(temp[2*i+1])

  # Make library
  totcites = 0
  for i in range(npapers):
    # If more than 7 authors the 8th author is  "..."
    # We now only extract the first 7 authors and if more "et al." is set as last author!
    authorlist = authors[i].split(",")
    if(len(authorlist)>= 8):
      authorlist[7] = "et al."
    googlecode = hrefs[i].split('citation_for_view=')[1]

    p = Paper(authors=authorlist,title=titles[i],journal=journal[i],cites=int(cites[i]),googlecode=googlecode)
    papers.append(p)
    totcites += p.cites

  # xxx We are still not adding abstract / all authors / papers that cite us... xxx
  # ...

  # Verbose we always provide
  print ''
  print '       Total papers in Google Scholar = ', npapers
  print '       Total cites  in Google Scholar = ', totcites
  print '       Cites per paper Google Scholar = ', float("{0:.2f}".format( totcites/float(npapers) ))

##############################################################################
##############################################################################

testinspire = True
testgoogle  = True
testnasa    = True
megaverbose = True

##############################################################################
##############################################################################
# Construct a library from an Inspire Hep id or name
##############################################################################
##############################################################################

if testinspire:
  inspire_hep_library = InspireHepLibrary(inspirename = "H.A.Winther.1")
  inspire_hep_library.construct_library(verbose = megaverbose)
  inspire_hep_library.compare_to_old_library("inspirelibrary.txt")
  inspire_hep_library.save_library("inspirelibrary.txt")
  inspire_hep_library.output_library_html("inspire_sample_output.html")

##############################################################################
##############################################################################
# Constuct a library from a Google Scholar id
##############################################################################
##############################################################################

if testgoogle:
  google_scholar_library = GoogleScholarLibrary(googleid = "oDCHqpAAAAAJ")
  google_scholar_library.construct_library(verbose = megaverbose)
  google_scholar_library.compare_to_old_library("googlelibrary.txt")
  google_scholar_library.save_library("googlelibrary.txt")
  google_scholar_library.output_library_html("google_sample_output.html")

##############################################################################
##############################################################################
# Construct a Nasa Ads library from a personal library or by searching
##############################################################################
##############################################################################

if testnasa:
  # nasa_ads_library = NasaAdsLibrary(get_by_search = True, searchphrase = "Winther, Hans A. ")
  nasa_ads_library = NasaAdsLibrary(get_personal_lib = True, nasaadslibid = ["My Papers", "5202f87fce"])
  nasa_ads_library.construct_library(verbose = megaverbose)
  nasa_ads_library.compare_to_old_library("library.txt")
  nasa_ads_library.save_library("library.txt")
  nasa_ads_library.output_library_html("sample_output.html")

