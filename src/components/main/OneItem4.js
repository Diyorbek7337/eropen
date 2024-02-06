import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';

const OneItem4 = () => {
    const { t } = useTranslation()
    const [ block, setBlock ] = useState(false)
    return (
        <div className="block__item">
                        <button className={block ? "block__title --active": "block__title"} type="button" data-spoller="" onClick={() => setBlock(!block)}>
                        <i className="fas fa-hand-holding-medical"></i> {t('Home.methodTitle')}

                          <span className="spoller__arrow"></span>
                        </button>
                        <div className={block ? "block__text --active": "block__text"} hidden="">
                          {t('Home.methodText')}
                           <br/>
                           {t('Home.methodText1')}
                        </div>
                      </div>
    );
};

export default OneItem4;