############################################################################
# Code to fetch information about scientific papers from arXiv, Nasa/ADS, 
# Google Scholar and/or Inspire HEP.
#
# * NASA/ADS: Fetch a personal ADS library or fetch all papers given searchwords
# * Inspire HEP: Fetch all papers given a HEP username / id
# * Google Scholar: Fetch all papers given a Google Scholar Profile id
# * arXiv: Fetch a list of recent papers in a category. Match against a list
#          of favoritte authors / keywords.
#
# Stores all papers we have fetched into a local library and output HTML with
# information about the papers. Stores the library to disc and reads it again
# every time we rerun the code and outputs informations about changes in the HTML.
#
# Code full of terrible handling of unicode/ascii/...!!
#
# Hans A. Winther (2015) (hans.a.winther@gmail.com)
##############################################################################

import time, sys, os.path, math                # Standard
import requests                                # Fetching webpages
import xml.etree.ElementTree as ET, lxml.etree # Parsing HTML/XML
import codecs

import itertools
from fuzzywuzzy import fuzz                    # String comparisons
 
# Useragent for arXiv fetches
ua_arxiv = {'Host':             'arxiv.org', 
            'User-Agent':       'Mozilla/5.0 (Macintosh; Intel Mac OS X 10.9; rv:39.0) Gecko/20100101 Firefox/39.0',
            'Accept':           'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language':  'en-GB,en;q=0.5',
            'Accept-Encoding':  'gzip, deflate'}

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

    # Temp variable for allowing reuse of methods...
    self.code = ""

##############################################################################
# Extracts the last 'n' new papers from arXiv in a given category e.g "astro-ph.CO"
# and returns a library with only title, arxivcode and authors
##############################################################################

def generate_url_arxiv(arxivcode, urltype = 'abs'):
  url = "http://arxiv.org/"+urltype+"/"+arxivcode
  return url

def extract_arXiv_recent(papers, category = "astro-ph", timeperiod = "pastweek", npapers = 100, verbose = False):
  url = "http://arxiv.org/list/"+category+"/"+timeperiod+"?skip=0&show="+str(npapers)
  response = requests.get(url,headers=ua_arxiv)
  html = response.text.encode('utf-8').strip()
  root = lxml.etree.HTML(html)

  arxivcodes = []
  for tag in root.findall(".//div[@id='content']/div[@id='dlpage']/dl/dt/span[@class='list-identifier']/a"):
    title = tag.get("title",None)
    if title == 'Abstract':
      arxivcodes.append(tag.text[6::].strip())

  titles = []
  for tag in root.findall(".//div[@id='content']/div[@id='dlpage']/dl/dd/div[@class='meta']/div[@class='list-title']"):
    titles.append(ET.tostring(tag,encoding='utf-8',method="text")[8::].strip())

  authors = []
  for tag in root.findall(".//div[@id='content']/div[@id='dlpage']/dl/dd/div[@class='meta']/div[@class='list-authors']"):
    tmp = []
    for a in tag.findall("a"):
      tmp.append(a.text)
    authors.append(tmp)

  na = len(authors)
  nt = len(titles)
  np = len(arxivcodes)

  print "We extracted",np,"new papers in arXiv:",category

  if(na != nt or na != np):
    print "Error we did not get all the data... (check read xml-format)",na,nt,np
    sys.exit(0)

  for i in range(np):
    p = Paper(title = titles[i], authors = authors[i], arxivcode = arxivcodes[i])
    papers.append(p)

  if verbose:
    for p in papers:
      print ''
      print '==================================================='
      print 'New recent paper on arXiv:',category
      print '==================================================='
      print p.arxivcode
      print p.title
      print " ; ".join(p.authors)

# String matching returns number in [0,1] according to how well the two strings are matched
def compare_strings(str1,str2):
  return fuzz.token_sort_ratio(str1, str2)/100.0

# Enter a list of names and search the library for authors. Outputs as html
def search_library_by_authors(filename, papers, myauthors, verbose = False):
  treshold = 0.85
  html = "" 
  indexlist = []

  for i,p in enumerate(papers):
    for a1 in p.authors:
      for a2 in myauthors:
        if compare_strings(a1,a2) > treshold: indexlist.append(i)
 
  indexlist = set(indexlist)
  for i in indexlist:
    if verbose:
      print ''
      print '==================================================='
      print 'Found new author in list...'  
      print '==================================================='
      print papers[i].arxivcode
      print papers[i].title
      print " ; ".join(papers[i].authors)
 
    # Fetch abstracts of papers...
    url = generate_url_arxiv(papers[i].arxivcode)
    response = requests.get(url,headers=ua_arxiv)
    xhtml = response.text.encode('utf-8').strip()
    root = lxml.etree.HTML(xhtml)
    for tag in root.findall(".//blockquote"):
      papers[i].abstract = ET.tostring(tag,encoding='utf-8',method="text").strip()

  # Make HTML with info about the papers we found
  html += "<h1>New papers by our favorite authors</h1>\n"
  html += "<p>\n"
  html += "  We searched for authors: " + " ; ".join(myauthors) + " and found " + str(len(indexlist)) + " new papers\n"
  html += "</p>\n"
  html += "<hr>"
  for i in indexlist:
    html += "<h2><a href=\""+ generate_url_arxiv(papers[i].arxivcode,'pdf') +"\">" + papers[i].title + "</a></h2>\n"
    html += "<p>\n"
    html += "  <b>Authors:</b> " + " ; ".join(papers[i].authors) + "\n"
    html += "</p>"
    html += "<p>\n"
    html += "  <b>Abstract:</b> " + papers[i].abstract + "\n"
    html += "</p>\n"
    html += "<hr>\n"

  with codecs.open(filename,'w',encoding='utf-8') as f:
    f.write(html)

##############################################################################
# NASA/ADS Library Class
##############################################################################

class NasaAdsLibrary:
  def __init__(self, get_personal_lib=False, nasaadslibid=None, get_by_search=False, searchphrase=None, name = "My Nasa/Ads Library"):
    self.get_personal_lib = get_personal_lib
    self.nasaadslibid     = nasaadslibid
    self.get_by_search    = get_by_search
    self.searchphrase     = searchphrase
    self.mylibrary        = []
    self.newciteinfo      = ""
    self.name             = name

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

  # Input: bibcode(s) for a paper. If more then one bibcode separate them by ';'
  def extract_bibtex(bibcodes, verbose = False):
    url  = generate_bib_url_NasaAds(bibcodes)
    html = get_html(url)
    allbibtex = html.split('@')[1::]
    bibtex = "@" + "@".join(allbibtex)
    return bibtex

##############################################################################
# Google Scholar Library Class
##############################################################################

class GoogleScholarLibrary:
  def __init__(self, googleid = None, name = "My Google Scholar Library"):
    self.googleid  = googleid
    self.mylibrary = []
    self.name      = name

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

def get_html_google(url):
  ua = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_8_2) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/27.0.1453.116 Safari/537.36'}
  response = requests.get(url, headers=ua)
  html = response.text.encode('utf-8').strip()
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

  with codecs.open(filename,'w',encoding='utf-8') as f:
    f.write(data)

##############################################################################
# Make HTML from a list of papers
# Input: a list of papers and the filename of the output
# Output: write filename to file and returns a string with the HTML code
##############################################################################

def whatisthis(s):
  if isinstance(s, str):
    return "string"
  elif isinstance(s, unicode):
    return "unicode"
  else:
    return "none"

def write_html_single_paper(paper,libtype,paperid):
  html  = "  <div class=\"paper\">\n"

  # Add title
  html += "    <h2><a href=\"#\" onclick=\"showhide('"+paperid+"')\">"+paper.title+"</a></h2>\n"
  html += "    <div class=\"paperinfo\" id=\"" + paperid + "\">\n"

  # Add authors (change to firstname first)
  html += "    <p>\n"
  html += "      <b>Authors:</b> "
  for author in paper.authors:
    if(libtype=="GoogleScholar"):
      if(whatisthis(author)=="unicode"):
        html += author.encode('ascii','ignore')
      else:
        html += author + " ; "
    else:
      html += author.split(',')[1] + " " + author.split(',')[0] + " ; "
  html += "\n"
  html += "    </p>\n"

  # Add abstract
  html += "    <p>\n"
  tmp = ""
  try:
    tmp += "      <b>Abstract:</b> "+paper.abstract+"\n"
  except:
    tmp  = "      <b>Abstract:</b> No abstract availiable\n"
    pass
  html += tmp
  html += "    </p>\n"

  # Add citation number
  html += "    <p>\n"
  html += "      <b>Citations (as of "+time.strftime("%d/%m/%Y")+"):</b> "+str(paper.cites)+"\n"
  html += "    </p>\n"

  # Add link to journal it is published in + arXiv + [Nasa ADS abstract / Inspire Hep / Google Scholar]
  html   += "    <p>\n    <b>Published in: </b>"
  if(paper.doi != ""):
    html += "<a href=\"http://dx.doi.org/"+paper.doi+"\">"+(paper.doi if paper.journal == "" else paper.journal)+"</a> | "
  if(paper.doi == "" and paper.journal != ""):
    html += paper.journal + " | "
  if(paper.arxivcode != ""):
    html += "<a href=\""+generate_url_arxiv(paper.arxivcode)+"\">eprint arXiv:"+paper.arxivcode+"</a> | "
  if(paper.bibcode != ""):
    html += "<a href=\""+generate_url_NasaAds(paper.bibcode,"ABSTRACT")+"\">Nasa ADS Abstract</a>\n"
  if(paper.inspirecode != ""):
    html += "<a href=\""+generate_url_InspireHep(paper.inspirecode, 'abstract')+"\">Inspire HEP Abstract</a>\n"
  if(paper.googlecode != ""):
    html += "<a href=\""+generate_url_GoogleScholar(paper.googlecode)+"\">Google Scholar Abstract</a>\n"
  html   += "    </p>\n"

  # Add links to paper that cite us
  if(paper.cites > 0):
    html += "    <p>\n"
    html += "      <b>Links to papers that cite us: </b>"
    for i,pp in enumerate(paper.cites_papers):
      if(libtype=="NasaAds"):
        html += "<a href=\"" + generate_url_NasaAds(pp.bibcode,"ABSTRACT") + "\">["+str(i+1)+"]</a>"
      elif(libtype=="InspireHep"):
        html += "<a href=\"" + generate_url_InspireHep(pp.inspirecode, 'abstract') + "\">["+str(i+1)+"]</a>"
      elif(libtype=="GoogleScholar"):
        html += "<a href=\"https://scholar.google.no/citations?view_op=view_citation&citation_for_view=" + pp.googlecode + "\">["+str(i+1)+"]</a>"
      html += " | " if i<paper.cites-1 else ""
    html += "\n"
    html += "    </p>\n"

  html += "    </div>\n"
  html += "    <hr>\n"
  html += "  </div>\n\n"

  return html

def write_library_as_html(papers, libname, newcites, filename, libtype, verbose):
  if(verbose):
    print ''
    print '==================================================='
    print 'Write information about library in HTML            '
    print '==================================================='

  # Sort papers by citation count
  papers.sort(key=lambda p: p.cites, reverse=True)

  # Use MathJax to process math in HTML?
  useMathJax = True

  # Calculate citation metrics
  totcites  = 0
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

  ############    START WRITE HTML     ###########
  html  = "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Draft//EN\">"
  html += "<html>\n"

  ############      WRITE HEADER       ###########
  html += "<head>\n\n"
  html += "  <title>"+libname+"</title>\n\n"
 
  # Add CSS
  html += "  <style type=\"text/css\">\n"
  html += "    .publications {display: hidden;}\n"
  html += "    .paper {display: hidden;}\n"
  html += "    .publicationinfo {display: hidden;}\n"
  html += "  </style>\n\n"

  # Add show/hide javascript
  html += "  <script type=\"text/javascript\">\n"
  html += "  function showhide(id) {\n"
  html += "    var e = document.getElementById(id);\n"
  html += "    if (e.style.display == 'block' || e.style.display==''){\n"
  html += "      e.style.display = 'none';\n"
  html += "    } else {\n"
  html += "     e.style.display = 'block';\n"
  html += "    }\n"
  html += "  }\n"
  html += "  </script>\n\n" 

  # Add MathJax Header
  if(useMathJax):
    html += "  <script type=\"text/x-mathjax-config\">\n"
    html += "    MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'],['\\\\(','\\\\)']]}});\n  </script>\n  <script type=\"text/javascript\" src=\"http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML\">\n"
    html += "  </script>\n\n"

  html += "</head>\n\n"
  ############    END WRITE HEADER     ###########

  ############    START WRITE BODY     ###########
  html += "<body>\n\n"

  # Show/hide papers link
  showhidelink="<a href=\"#\" onclick=\""
  for i, paper in enumerate(papers):
    paperid = "paper_" + str(i+1)
    showhidelink +="showhide('"+paperid+"');"
  showhidelink +="\">Show/hide all papers</a>"

  # Add introduction with some useful metrics
  html += "<div class=\"libraryinfo\">\n"
  html += "  <h1>"+libname+"</h1>\n"
  html += "  <p>\n"
  html += "    We found a total of "+str(totpapers)+" papers with a total of "+str(totcites)+" citations. This corresponds to "+str(citeperpaper)+" citations per paper. The total number of collaborators on these papers is " +str(ncollaborator)+".\n    <br>\n"
  html += "    " + showhidelink + " | " + "<a href=\"#\" onclick=\"showhide('publicationinfo')\">Show/hide newinfo</a>\n"
  html += "  </p>\n"
  html += "  <hr>\n"
  html += "</div>\n\n"

  ############     WRITE NEWINFO      ###########
  # Add info about changes since last time we ran this code
  html += "<div class=\"publicationinfo\" id=\"publicationinfo\">\n"
  html += "  <h2>New papers and citations</h2>\n"
  if(newcites != ""):
    html += newcites
  else:
    html += "  <p>\n    No new papers/citations since last time we checked.\n  </p>\n"
  html += "  <hr>\n"
  html += "</div>\n\n"
  ############   END WRITE NEWINFO    ###########

  ############   WRITE PUBLICATIONS   ###########

  # Loop over all papers in library and write info
  html += "<div class=\"publications\" id=\"publications\">\n"
  for i, paper in enumerate(papers):
    paperid = "paper_" + str(i+1)
    html += write_html_single_paper(paper,libtype,paperid)
  html += "</div>\n\n"

  ############ END WRITE PUBLICATIONS ###########

  html += "</body>\n"
  html += "</html>"

  ############     END WRITE HTML      ###########

  # Encode to ascii using xml characters
  if(libtype=="NasaAds"):
    html = html.encode('ascii', 'xmlcharrefreplace')
   
  # xxx Decode to unicode
  html = html.decode('utf-8')

  if(verbose):
    print "Write output to : ", filename
  with codecs.open(filename,'w',encoding='utf-8') as f:
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
  authors = sorted(set(authors)) # Not needed
  nauthors = len(authors)

  for a,b in itertools.combinations(authors,2):
    if(compare_strings(a,b) > 0.75):
      print "Author",a,"is equal to",b,"corr=",compare_strings(a,b)
      if(a in authors): authors.remove(a)
  return authors

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
  url = "http://inspirehep.net/author/profile/citations-summary"
  form_data = {'jsondata': '{\"personId\":'+str(inspireid)+'}',}
  html = get_html(url, form_data)
  root = lxml.etree.HTML(html)

  # Citations
  tags      = root.findall(".//tr/td")
  cites     = int(tags[3].text)
  cites_pub = int(tags[4].text)

  # Number of papers
  tags        = root.findall(".//table/tbody/tr/td/a")  
  npapers     = int(tags[0].text)
  npapers_pub = int(tags[1].text)
 
  # Cite per paper
  citeperpaper     = 0.0 if npapers     == 0 else float("{0:.2f}".format( cites/(float(npapers)) ))
  citeperpaper_pub = 0.0 if npapers_pub == 0 else float("{0:.2f}".format( cites_pub/(float(npapers_pub)) ))

  # Fetch papers (Inspire only display 250 at the time)
  for i in range( int( math.ceil( abs(npapers-0.5)/250.0 ) ) ):
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

  print ''
  print '       Total papers in Inspire HEP    = ', npapers,'(',npapers_pub,')'
  print '       Total cites  in Inspire HEP    = ', cites,'(',cites_pub,')'
  print '       Cites per paper Inspire HEP    = ', citeperpaper,'(',citeperpaper_pub,')' 

##############################################################################
# Extracts title, abstract, authors and doi (if it exists) from a inspirecode
# and stores it in paper
##############################################################################

def extract_main_paper_data_InspireHep(paper, inspirecode, verbose):
    url = generate_url_InspireHep(paper.inspirecode,"xml")
    xml = get_html(url)
    root = lxml.etree.HTML(xml) 
 
    # xxx should not be needed but it is for some reason!
    paper.authors = []
    paper.title = ""
    paper.abstract = ""

    if(verbose):
      print ''
      print 'Extracting data for inspirecode = ', inspirecode
      print 'Url = ', url

    for author   in root.iter("author"):
      if author.text   != None: paper.authors.append(author.text.encode('ascii', 'xmlcharrefreplace'))
      if author.text != None: print author.text.encode('ascii', 'xmlcharrefreplace')
    for abstract in root.iter("abstract"):
      if abstract.text != None: paper.abstract = abstract.text.encode('ascii', 'xmlcharrefreplace')
    for title    in root.iter("title"):
      if title.text    != None: paper.title = title.text.encode('ascii','xmlcharrefreplace')
    for doi      in root.iter("electronic-resource-num"):
      if doi.text      != None: paper.doi = doi.text

    if(verbose):
      print paper.title,'|',paper.authors,'|',paper.doi

##############################################################################
# Generate different url's related to an Inspire HEP entry
##############################################################################

def generate_url_InspireHep(inspirecode, urltype = ''):
  # Map between code and url postfix
  dictionary = {'bibtex': 'export/hx', 'latex': 'export/hlxu', 'xml': 'export/xe', 'harvmac': 'export/hlxh', 'abstract': '', '': '', 'citations': 'citations', 'marcxml': 'export/xm', 'marc': 'export/hm', 'nml': 'export/xn','dc': 'export/xd'}
  code = urltype
  for key in dictionary:
    code = code.replace(key,dictionary[key])
  url = "http://inspirehep.net/record/"+inspirecode+"/"+code
  return url

##############################################################################
# Collect information about all papers that cite a given paper
##############################################################################

def get_inspirecodes_of_cites(inspirecode, verbose):
  # Init
  cites_inspirecode = []; cites = 0

  if(verbose):
    print ''
    print '==================================================='
    print 'Fetching inspirecodes of papers that cite',inspirecode
    print '==================================================='

  # Loop over different pages containing data (since 250 max per page)
  for i in range(100):
    url  = "http://inspirehep.net/search?ln=en&ln=en&p=refersto:recid:"+inspirecode+"&of=hb&action_search=Search&rg=250&jrec="+str(i*250+1)
    html = get_html(url)
    root = lxml.etree.HTML(html)

    # Get number of cites
    if i==0:
      for c in root.findall(".//td[@class='searchresultsboxheader'][@align='center']/strong"):
        cites = int(c.text)

    # Get all inspirecodes from this page
    for a in root.findall(".//div[@class='record_body']/a[@class='titlelink']"):
      cites_inspirecode.append(a.get("href",None).split("record/")[1])
 
    # Exit when we have fetched all papers xxx NB: Check this. Should it really be cites here and not npapers?
    if(i >= int( math.floor( abs(cites-0.5)/250.0 ) ) ): break

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
  oldpapers = []; html = ""

  # Copy to paper-code to the temp 'code' variable
  for p in papers:
    p.code = papercode_from_libtype(p, libtype)
    for c in p.cites_papers:
      c.code = papercode_from_libtype(c, libtype)
  
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
  with codecs.open(filename,'r',encoding='utf-8') as f:
    data = f.readlines()
    if data == []: return html

    noldpapers = int(data[0])
    for i in range(noldpapers):
      code = data[2+5*i].strip()
      title = data[3+5*i].strip()
      cites = int(data[4+5*i])
      cite_codes = data[5+5*i].strip().split(';')
      q = Paper(title=title,cites=cites)
      q.code = code

      # Add papers that cite us (only code)
      for cb in cite_codes:
        newpaper = Paper()
        newpaper.code = cb
        q.cites_papers.append(newpaper)
      oldpapers.append(q)

  # How many new papers?
  newpapers = len(papers) - len(oldpapers)
  if(verbose): print 'Comparing old library with current one. We have',newpapers,'new papers'
  if(newpapers>0):
    html += "<p>\n"
    html += "  We have "+str(newpapers)+" papers\n"
    html += "</p>\n"

  # Loop over all papers and compare old and new library [only check papers and codes of cites]
  newcites = 0
  for newp in papers:
    allnewcitecode = []

    # Check if paper newp is also in old library
    match = [oldp for oldp in oldpapers if oldp.code == newp.code]

    if(len(match)==0):
      # Paper newp is not in oldpapers
      if(verbose):
        print "We found a new paper:", newp.code, "that is not in the old library"
        print "This paper has",newp.cites,"citations"
      newcites += newp.cites
      for p in newp.cites_papers:
        allnewcitecode.append(p.code)
    else:
      # Paper newp is in oldpapers. Check if citations has changed
      oldp = match[0]
      if(newp.cites > oldp.cites):
        if(verbose): print "Paper",newp.code,"have",newp.cites-oldp.cites,"new citations!"
        newcites += newp.cites-oldp.cites

        # Find bibcodes of the new citations
        oldcite_code = []; newcite_code = []
        for p in oldp.cites_papers:
          oldcite_code.append(p.code)
        for p in newp.cites_papers:
          newcite_code.append(p.code)
        for code in newcite_code:
          if code not in oldcite_code:
            allnewcitecode.append(code)

    # Compile some info about the new cites in HTML
    if(allnewcitecode != []):
      html += "<p>\n"
      html += "  New citations for paper <b>\""+newp.title+"\"</b>: "
      for code in allnewcitecode:
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
  form_data = {'bibcode': bibcodes, 'sort': 'CITATIONS', 'data_type': 'XML', 'submit': 'submit',}
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
  ns = {'ref': 'http://ads.harvard.edu/schema/abs/1.1/abstracts'}
  if(verbose):
    print ''
    print '==================================================='
    print 'Adding all papers to library                       '
    print '==================================================='

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

    # Add paper to library
    p = Paper(title=title,authors=authors,bibcode=bibcode,cites=cites,abstract=abstract)
    extract_publication_info_NasaAds(p, verbose)
    papers.append(p)
    if(verbose): print '  * Adding paper: \"',p.title,"\" Cites: ", p.cites

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
  bibcodes = []
  url = generate_url_NasaAds(bibcode, "CITATIONS")
  if(verbose):
    print ''
    print '==================================================='
    print 'Get bibcodes of papers that cite', bibcode
    print '==================================================='
    print "Url = ", url

  html = get_html(url)
  root = lxml.etree.HTML(html)
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
  doi  = arxiv = None

  # Fetch meta tags from abtract page
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
      # Get all info about paper that cite us
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
# Make Google Scholar requests. As done now we only fetch the first 100 papers 
# on the page! Addition requests (form submits) are needed to extract the rest!
##############################################################################

def construct_library_GoogleScholar(papers, googleid, verbose):
  url = "https://scholar.google.com/citations?user="+googleid
  hrefs = []; authors = []; journal = []; cites = []; temp = []; titles = []
  
  if(verbose):
    print ''
    print '==================================================='
    print 'Construct Google Scholar library', googleid
    print '==================================================='
    print 'Fetch url:', url

  # Extract data from HTML
  html = get_html_google(url)
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

  # Make library...
  totcites = 0
  for i in range(npapers):
    googlecode = hrefs[i].split('citation_for_view=')[1]
    p = Paper(title=titles[i],journal=journal[i],cites=int(cites[i]),googlecode=googlecode)
    papers.append(p)
    extract_main_paper_data_GoogleScholar(p, googlecode, verbose)
    totcites += p.cites

  # NB: We are still not adding papers that cite us!
  # ...

  print ''
  print '       Total papers in Google Scholar = ', npapers
  print '       Total cites  in Google Scholar = ', totcites
  print '       Cites per paper Google Scholar = ', float("{0:.2f}".format( totcites/float(npapers) ))

##############################################################################
# Extracts journal/abstract/cites from a GoogleScholar abstract page
##############################################################################

def generate_url_GoogleScholar(googlecode):
  url = "https://scholar.google.no/citations?view_op=view_citation&hl=en&citation_for_view="+googlecode
  return url

def lookup_entry_in_list_GoogleScholar(mylist, entry, errorvalue = ""):
  value = filter(lambda field: field[0] == entry, mylist)
  return value[0][1] if value != [] else errorvalue

def extract_main_paper_data_GoogleScholar(paper, googlecode, verbose):
  url  = generate_url_GoogleScholar(googlecode)
  html = get_html_google(url)
  root = lxml.etree.HTML(html)

  if(verbose):
    print 'Extracting data for paper',googlecode
    print 'Url = ',url

  # Fetch data from a Google Scholar abstract page
  temp = []
  for tag in root.findall(".//div[@class='gs_scl']"):
    field = None; value = None
    # Get field description
    for tag2 in tag.findall(".//div[@class='gsc_field']"):
      field = tag2.text
    # Get field value
    if(field == "Description"):
      # For the abstract field
      for tag2 in tag.findall(".//div[@id='gsc_descr'][@class='gsc_value']"):
        value = ET.tostring(tag2,method="text",encoding="utf-8")[9::] # Extract xx in "Abstract xx"
    elif(field == "Total citations"):
      for tag2 in tag.findall(".//div[@class='gsc_value']/div/a"):
        value = tag2.text[8::] # Extract xx in "Cited by xx"
    else:
      for tag2 in tag.findall(".//div[@class='gsc_value']"):
        value = tag2.text
    # Add info to list
    if (field != None and value != None):
      temp.append([field, value])

  # Extract information
  cites    = int(lookup_entry_in_list_GoogleScholar(temp, "Total citations", "0"))
  abstract =     lookup_entry_in_list_GoogleScholar(temp, "Description")
  journal  =     lookup_entry_in_list_GoogleScholar(temp, "Journal", "[Journal not provided]")
  volume   =     lookup_entry_in_list_GoogleScholar(temp, "Volume")
  pages    =     lookup_entry_in_list_GoogleScholar(temp, "Pages")
  issue    =     lookup_entry_in_list_GoogleScholar(temp, "Issue")
  authors  =     lookup_entry_in_list_GoogleScholar(temp, "Authors").split(',')
  journal  = (journal+" "+volume+" "+issue+" "+pages).strip()
  for i,author in enumerate(authors): authors[i] = author.strip()

  # Add information to paper
  paper.cites    = cites
  paper.journal  = journal
  paper.authors  = authors
  paper.abstract = abstract

  if(verbose):
    print '* Adding paper data: Cites=',cites,"Journal=",journal,"Authors=",";".join(authors)

##############################################################################
##############################################################################

runarxiv    = True
runinspire  = True
rungoogle   = True
runnasa     = True
megaverbose = True

##############################################################################
##############################################################################
# Construct a library from an Inspire Hep id or name
##############################################################################
##############################################################################

if runinspire:
  inspire_hep_library = InspireHepLibrary(inspirename = "H.A.Winther.1")#"S.Kirkevold.Naess.1") #"H.A.Winther.1")
  inspire_hep_library.construct_library(verbose = megaverbose)
  inspire_hep_library.compare_to_old_library("inspirelibrary.txt")
  inspire_hep_library.save_library("inspirelibrary.txt")
  inspire_hep_library.output_library_html("inspire.html")

##############################################################################
##############################################################################
# Constuct a library from a Google Scholar id
# NB: Google can block requests if one ties too many in a short enough time
# This gives Error: 503
##############################################################################
##############################################################################

if rungoogle:
  google_scholar_library = GoogleScholarLibrary(googleid = "oDCHqpAAAAAJ")
  google_scholar_library.construct_library(verbose = megaverbose)
  google_scholar_library.compare_to_old_library("googlelibrary.txt")
  google_scholar_library.save_library("googlelibrary.txt")
  google_scholar_library.output_library_html("google.html")

##############################################################################
##############################################################################
# Construct a Nasa Ads library from a personal library or by searching
##############################################################################
##############################################################################

if runnasa:
  # nasa_ads_library = NasaAdsLibrary(get_by_search = True, searchphrase = "Winther, Hans A. ")
  nasa_ads_library = NasaAdsLibrary(get_personal_lib = True, nasaadslibid = ["My Papers", "5202f87fce"])
  nasa_ads_library.construct_library(verbose = megaverbose)
  nasa_ads_library.compare_to_old_library("nasaadslibrary.txt")
  nasa_ads_library.save_library("nasaadslibrary.txt")
  nasa_ads_library.output_library_html("nasaads.html")

##############################################################################
##############################################################################
# Fetch recent arxiv papers and search for authors and output to HTML
# what we find
##############################################################################
##############################################################################

if runarxiv:
  arxivpapers = []
  extract_arXiv_recent(arxivpapers, "astro-ph", "pastweek", 1000, megaverbose)
  myauthors = ["Hans A. Winther", "Sigurd Naess", "Philippe Brax", "David F. Mota", "Pedro Ferreira", "David Alonso", "Baojiu Li", "Fabian Schmidt", "Luca Amendola", "Kazuya Koyama", "Alexandre Barreira"]
  search_library_by_authors("arxiv.html", arxivpapers, myauthors, megaverbose)

