def getCondition(region='Vrancea'):

    # VRANCEA
    if region=='Vrancea':
        year = input('Input year from which you want visualization (recommended: 1976) : ')
        dateandtime = year + '-01-01 00:00:00'
        condition=(f"SELECT * FROM romplus WHERE `dateandtime`>='{dateandtime}'"
            f" AND `latitude`>=45.2 AND `latitude`<=46 AND `longitude`>=26" 
            f" AND `longitude`<=27 AND `depth`>=50 AND `depth`<=200")
        
    # ROMANIA
    if region=='Romania':
        year = input('Input year from which you want visualization (recommended: 1976) : ')
        dateandtime = year + '-01-01 00:00:00'
        condition=f"SELECT * FROM romplus WHERE `dateandtime`>='{dateandtime}'"

    # CALIFORNIA
    if region=='California':
        year = input('Input year from which you want visualization (recommended: 1984) : ')
        dateandtime = year + '-01-01 00:00:00'
        condition=(f"SELECT * FROM california WHERE `dateandtime`>='{dateandtime}'"
                    f" AND (`magtype` LIKE 'l' OR `magtype` LIKE 'w')")

    # ITALY
    if region=='Italy':
        year = input('Input year from which you want visualization (recommended: 1986) : ')
        dateandtime = year + '-01-01 00:00:00'
        condition=f"SELECT * FROM italy WHERE `dateandtime`>='{dateandtime}' AND `depth`>=0"

    # JAPAN
    if region=='Japan':
        year = input('Input year from which you want visualization (recommended: 1992) : ')
        dateandtime = year + '-01-01 00:00:00'
        condition=(f"SELECT * FROM japan WHERE `dateandtime`>='{dateandtime}'"
                    f" AND `latitude`>0 AND `longitude`>0 AND `depth`>0"
                    f" AND `magtype` LIKE 'V'")

    return condition,year