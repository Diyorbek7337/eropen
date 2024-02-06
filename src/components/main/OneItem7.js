import React, { useState } from 'react';
import { useTranslation } from 'react-i18next';

const OneItem7 = () => {
    const { t } = useTranslation()
    const [ block, setBlock ] = useState(false)
    return (
        <div className="block__item">
                        <button className={block ? "block__title --active": "block__title"} type="button" data-spoller="" onClick={() => setBlock(!block)}>
                        <i className="fas fa-business-time"></i>{t('Home.storageTitle')}
                          <span className="spoller__arrow"></span>
                        </button>
                        <div className={block ? "block__text --active": "block__text"} hidden="">
                        {t('Home.storageText1')}
                          <br/>
                          {t('Home.storageText2')}<br/>
                        </div>
                      </div>
    );
};

export default OneItem7;