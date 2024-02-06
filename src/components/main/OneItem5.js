import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';

const OneItem5 = () => {
    const { t } = useTranslation()
    const [ block, setBlock ] = useState(false)
    return (
        <div className="block__item">
                        <button className={block ? "block__title --active": "block__title"} type="button" data-spoller="" onClick={() => setBlock(!block)}>
                        <i className="fa fa-comment-medical"></i> {t('Home.sideTitle')}
                          
                          <span className="spoller__arrow"></span>
                        </button>
                        <div className={block ? "block__text --active": "block__text"} hidden="">
                           {t('Home.sideText1')}
                           <br/>
                           {t('Home.sideText2')} <br/>
                           {t('Home.sideText3')}
                        </div>
                      </div>
    );
};

export default OneItem5;