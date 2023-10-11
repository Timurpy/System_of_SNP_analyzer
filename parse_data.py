from bs4 import BeautifulSoup as BS
import requests

def get_MAF_info_from_dbSNP(rs_id, MAF_threshold):
    url = f"https://www.ncbi.nlm.nih.gov/snp/{rs_id}"
    # Page request and HTML analysis using BeautifulSoup
    with requests.Session() as s:
        response = s.get(url)
        soup = BS(response.text, "html.parser")
        
        try: 
            info = soup.find_all("dd")[4]
        except IndexError:
            return []
        
        info = str(info).split('\n')
        info = [i.lstrip().replace('</div>', '').replace('<div>', '').replace('<div>', '') for i in info if (not '<span>' in i)]
        info = [i for i in info if not i == '']
        allele_data = []
        i = 1
        while i < len(info) - 1:
            if '(' in info[i] and ')' in info[i]:
                allele_data.append(info[i])
                i += 1
            else:
                allele_data.append(f"{info[i]} {info[i+1]}")
                i += 2
        allele_data = sorted(allele_data, key=lambda x: float(x.split()[0].split('=')[-1]), reverse=True)
        allele_data = [i for i in allele_data if float(i.split()[0].split('=')[-1]) >= MAF_threshold]
        return allele_data