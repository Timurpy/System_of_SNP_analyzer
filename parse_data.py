from bs4 import BeautifulSoup as BS
import requests

def get_info_from_dbSNP(rs_id, MAF_threshold):
    url = f"https://www.ncbi.nlm.nih.gov/snp/{rs_id}"
    # Page request and HTML analysis using BeautifulSoup
    with requests.Session() as s:
        response = s.get(url)
        soup = BS(response.text, "html.parser")
        clinical_info = soup.find("table", {'id': 'clinical_significance_datatable'})
        
        # scrapin maf info
        try: 
            maf_info = soup.find_all("dd")[4]
        except IndexError:
            return [], None
        
        maf_info = str(maf_info).split('\n')
        maf_info = [i.lstrip().replace('</div>', '').replace('<div>', '').replace('<div>', '') for i in maf_info if (not '<span>' in i)]
        maf_info = [i for i in maf_info if not i == '']
        allele_data = []
        i = 1
        while i < len(maf_info) - 1:
            if '(' in maf_info[i] and ')' in maf_info[i]:
                allele_data.append(maf_info[i])
                i += 1
            else:
                allele_data.append(f"{maf_info[i]} {maf_info[i+1]}")
                i += 2
        
        # scrapin clinical significance info
        try:
            clinical_info = clinical_info.find("tbody")
            clinical_info = clinical_info.find_all("tr")
            
            clinical_reports = []
            for entry in clinical_info:
                entry_info = entry.find_all("td")
                disease_name = entry_info[-2].text
                clinical_significance = entry_info[-1].text
                clinical_report =f"{disease_name}: {clinical_significance}"
                clinical_reports.append(clinical_report)
            final_clinical_info = '(' + '; '.join(clinical_reports) + ')'

        except AttributeError:
            final_clinical_info = None
            pass
        
        
        allele_data = sorted(allele_data, key=lambda x: float(x.split()[0].split('=')[-1]), reverse=True)
        allele_data = [i for i in allele_data if float(i.split()[0].split('=')[-1]) >= MAF_threshold]
        return allele_data, final_clinical_info