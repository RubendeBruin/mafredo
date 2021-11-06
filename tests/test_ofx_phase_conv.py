from mafredo import Hyddb1
import matplotlib.pyplot as plt

def test_read_ofx_conventions():
    lags = Hyddb1.create_from_orcaflex_yml(filename=r'files/refcase_ofx_lags.yml', vessel_type_name='Vessel type1')
    leads = Hyddb1.create_from_orcaflex_yml(filename=r'files/refcase_ofx_leads.yml', vessel_type_name='Vessel type1')
    leads_zd = Hyddb1.create_from_orcaflex_yml(filename=r'files/refcase_ofx_leads_zero_down.yml', vessel_type_name='Vessel type1')

    rao90_lags = lags.force_rao(3) # roll
    rao90_leads = leads.force_rao(3)
    rao90_leads_zd = leads_zd.force_rao(3)

    rao90_lags.regrid_direction(90)
    rao90_leads.regrid_direction(90)
    rao90_leads_zd.regrid_direction(90)

    rao90_lags._data['phase'].plot(marker='*')
    rao90_leads._data['phase'].plot(marker='o')
    rao90_leads_zd._data['phase'].plot()

    plt.show()
