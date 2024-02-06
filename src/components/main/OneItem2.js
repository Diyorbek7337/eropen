
import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';
const OneItem2 = () => {
    const { t } = useTranslation()
    const [ block, setBlock ] = useState(false)
    return (
        <div className="block__item">
                        <button className={block ? "block__title --active": "block__title"} type="button" data-spoller="" onClick={() => setBlock(!block)}>
                        <i className="fa fa-capsules"></i> {t('Home.pharmacologicalTitle')}

                          <span className="spoller__arrow"></span>
                        </button>
                        <div className={block ? "block__text --active": "block__text"} hidden="">
                          {t('Home.pharmacologicalText1')}
                          <br/>
                          {t('Home.pharmacologicalText2')}
                           <br/>
                           {t('Home.pharmacologicalText3')}
                           <br/>
                           {t('Home.pharmacologicalText4')}
                           <br/>
                          {t('Home.pharmacologicalText5')}
                          <br/>
                         {t('Home.pharmacologicalText6')}
                         <br/>
                        </div>
                      </div>
    );
};

export default OneItem2;